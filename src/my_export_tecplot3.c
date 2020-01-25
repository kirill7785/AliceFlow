// Файл my_export_tecplot3.c передача результатов
// моделирования в программу tecplot360

#pragma once
#ifndef MY_EXPORT_TECPLOT3_C
#define MY_EXPORT_TECPLOT3_C 1

//#include "windows.h" // для функции WinExec
#include <string.h>


// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
void exporttecplotxy360_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal** potent, doublereal **rhie_chow)
{
	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	fp = fopen64("ALICEFLOW0_03.PLT", "w");
#else
	err = fopen_s(&fp, "ALICEFLOW0_03.PLT", "w");
#endif

	
	// создание файла для записи.
	if ((err) != 0) {
		printf("Create File Error\n");
		exit(1);
	}
	else {
		if (fp != NULL) {
			// запись заголовка
			fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

			// запись имён переменных
			fprintf(fp, "VARIABLES = x, y, z, Vx, Vy, Vz, Mag, Pressure , Normal , PAm, Zero\n");

			// запись информации о зонах
#if doubleintprecision == 1
			//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
			//if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
			fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
			//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
			//if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
			fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif
			

			integer i = 0; // счётчики 
			//integer j = 0; // цикла for

			// запись x
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись y
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись z
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");


			// запись горизонтальной Vx скорости
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[VX][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись вертикальной Vy скорости
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[VY][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			fprintf(fp, "\n");

			// запись Vz скорости
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[VZ][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			fprintf(fp, "\n");


			// запись модуля скорости
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", sqrt(potent[VX][i] * potent[VX][i] + potent[VY][i] * potent[VY][i] + potent[VZ][i] * potent[VZ][i]));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись давления
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[PRESS][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись Normal
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", rhie_chow[0][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись PAm
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[PAM][i]); // rhie_chow[1][i]
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись Rhi-Chow
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", rhie_chow[2][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись информации о разностной сетке
			for (i = 0; i < ncell; i++) {
#if doubleintprecision == 1
				if (ionly_solid_visible == 0) {
					fprintf(fp, "%lld %lld %lld %lld ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			}
				if (ionly_solid_visible == 1) {

					fprintf(fp, "%lld %lld %lld %lld ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
				}
#else
				if (ionly_solid_visible == 0) {
					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
				}
				if (ionly_solid_visible == 1) {

					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
				}
#endif
				
			}

			if (fp != nullptr) {
				fclose(fp); // закрытие файла
			}
			else {
				printf("exporttecplotxy360_3D problem close file .\n");
				system("pause");
				exit(1);
			}
			printf("file is successfully written...OK.\n");
			// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
			   //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL
		}
	}
} // exporttecplotxy360_3D

// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
void exporttecplotxy360T_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal* potent)
{
	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	fp = fopen64("ALICEFLOW0_03.PLT", "wb");
#else
	err = fopen_s(&fp, "ALICEFLOW0_03.PLT", "wb");
#endif
	
	// создание файла для записи.
	if ((err) != 0) {
		printf("Create File Error\n");
	}
	else {
		if (fp!=NULL) {

		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

		// запись имён переменных
		fprintf(fp, "VARIABLES = x, y, z, Temp\n");

#if doubleintprecision == 1
		// запись информации о зонах
		fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
		// запись информации о зонах
		fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif

		

		integer i = 0; // счётчики 
		//integer j = 0; // цикла for

		// запись x
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись y
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись z
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");


		// запись температуры
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", potent[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись информации о разностной сетке
		for (i = 0; i < ncell; i++) {
#if doubleintprecision == 1
			//fprintf(fp, "%lld %lld %lld %lld ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
			//fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i], nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#else
			//fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
			//fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			fprintf(fp, "%d %d %d %d %d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i], nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#endif

				}

		fclose(fp); // закрытие файла
		printf("file is successfully written...OK.\n");
		// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
		   //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);
	}
	}
} // exporttecplotxy360T_3D

// Трёхэтапная запись файла позволяет сэкономить 19N оперативной памяти.
// Далее везде применяется трёхэтампное формирование выходного файла.

// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// части 1 и 3.
void exporttecplotxy360T_3D_part1and3(TEMPER &t, integer maxelm, integer maxbound, bool bextendedprint, integer ncell,
									  integer** nvtx, integer** nvtxcell, TOCHKA* pa,
									  BOUND* sosedb, integer ivarexport, integer** ptr_out)
{

	if (bvery_big_memory) {

		t.database.maxelm = maxelm;
		t.database.ncell = ncell;

		// extended printeger не предусмотрено.

		// Если память уже выделялась ранее то её надо освободить.
		if (t.database.x != nullptr) {
			free(t.database.x);
			t.database.x = nullptr;
		}
		if (t.database.y != nullptr) {
			free(t.database.y);
			t.database.y = nullptr;
		}
		if (t.database.z != nullptr) {
			free(t.database.z);
			t.database.z = nullptr;
		}
		if (t.database.ptr != nullptr) {
			for (integer i_92 = 0; i_92 < 2; i_92++) {
				delete[] t.database.ptr[i_92];
				t.database.ptr[i_92] = nullptr;
			}
			delete[] t.database.ptr;
			t.database.ptr = nullptr;
		}
		if (t.database.nvtxcell != nullptr) {
			for (integer i_92 = 0; i_92 < 8; i_92++) {
				delete[] t.database.nvtxcell[i_92];
				t.database.nvtxcell[i_92] = nullptr;
			}
			delete[] t.database.nvtxcell;
			t.database.nvtxcell = nullptr;
		}

		t.database.x = (doublereal*)malloc(maxelm*sizeof(doublereal));
		t.database.y = (doublereal*)malloc(maxelm*sizeof(doublereal));
		t.database.z = (doublereal*)malloc(maxelm*sizeof(doublereal));
		t.database.ptr = new integer*[2];
		for (integer j = 0; j < 2; j++) {
			t.database.ptr[j] = new integer[maxelm];
		}
		for (integer j = 0; j < 2; j++) {
			for (integer i = 0; i < maxelm; i++) {
				t.database.ptr[j][i] = ptr_out[j][i];
			}
		}
		t.database.nvtxcell = new integer*[8];
		for (integer i = 0; i < 8; i++) {
			t.database.nvtxcell[i] = new integer[ncell];
		}

		// запись x
		if (t.database.x != nullptr) {
			for (integer i = 0; i < maxelm; i++) {
				t.database.x[i] = 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x);
			}
		}
		else {
			printf("ERROR malloc: t.database.x ==nullptr. exporttecplotxy360T_3D_part1and3 \n");
			system("pause");
			exit(1);
		}
		// запись y
		if (t.database.y != nullptr) {
		    for (integer i = 0; i < maxelm; i++) {
			    t.database.y[i] = 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y);
		    }
	    }
	    else {
	    	printf("ERROR malloc: t.database.y ==nullptr. exporttecplotxy360T_3D_part1and3 \n");
	    	system("pause");
	     	exit(1);
	    }

		// запись z
		if (t.database.z != nullptr) {
			for (integer i = 0; i < maxelm; i++) {
				t.database.z[i] = 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z);
			}
		}
		else {
			printf("ERROR malloc: t.database.z ==nullptr. exporttecplotxy360T_3D_part1and3 \n");
			system("pause");
			exit(1);
		}

		// запись информации о разностной сетке
		for (integer i = 0; i < ncell; i++) {
			for (integer j = 0; j < 8; j++) {
				t.database.nvtxcell[j][i] = nvtxcell[j][i];
			}
		}

	}
	else {

		t.database.maxelm = 0;
		t.database.ncell = 0;
		//Запись в файл очень медленная.
		t.database.x = nullptr;
		t.database.y = nullptr;
		t.database.z = nullptr;
		t.database.nvtxcell = nullptr;
		t.database.ptr = nullptr;


		// расширенная печать
		// При расширенной печати мы печатаем также и граничные узлы.
		// bextendedprint=true; расширенная печать. 

		// ivarexport == 1 печатается только поле температур,
		// ivarexport == 2 печатается только гидродинамика,
		// ivarexport == 3 печатается и поле температур и гидродинамика.

		FILE *fp=NULL;
		errno_t err=0;
#ifdef MINGW_COMPILLER
		fp = fopen64("ALICEFLOW0_06_temp_part1.txt", "w");
#else
		err = fopen_s(&fp, "ALICEFLOW0_06_temp_part1.txt", "w");
#endif
		
		// создание файла для записи:
		// файл состоит из трёх частей: 
		// 1 и 3 часть записываются сразу
		// вторая часть с результатами расчёта записывается
		// после расчёта. Такая трёхэтапная запись файла выбрана в целях
		// сокращения объёма используемой оперативной памяти.
		// Экономия памяти 19N.

		// запись частей 1 и 3
		if ((err) != 0) {
			printf("Create File temp part1 Error\n");
			//getchar();
			system("pause");

		}
		else {
			// запись заголовка
			/*fprintf(fp, "TITLE = \"ALICEFLOW0_06\"\n");

			// запись имён переменных
			switch (ivarexport) {
			case 1 : fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n"); break; // печатается только поле температур
			case 2 : fprintf(fp, "VARIABLES = x, y, z, Speed, Pressure, Vx, Vy, Vz, Rho, Mu, Mut\n"); break; // печатается только гидродинамика
			case 3 : fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, Speed, Pressure, Vx, Vy, Vz, Rho, Mu, Mut\n"); break; // печатается и температура и гидродинамика
			default : printf("Error export tecplot. Nonselected exporting variables...\n"); getchar();
			}

			#if doubleintprecision == 1
			    // запись информации о зонах
			    fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
			#else
			    // запись информации о зонах
			    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
			#endif

			
			*/

			integer i = 0; // счётчики 
			//integer j=0; // цикла for

			// запись x
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (sosedb[i].Norm) {// определим внутреннюю нормаль к границе
					case ESIDE: fprintf(fp, "%+.16f ", pa[nvtx[0][sosedb[i].iI] - 1].x);
						break;
					case WSIDE: fprintf(fp, "%+.16f ", pa[nvtx[1][sosedb[i].iI] - 1].x);
						break;
					case NSIDE : fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].x + pa[nvtx[1][sosedb[i].iI] - 1].x));
						break;
					case SSIDE:fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].x + pa[nvtx[1][sosedb[i].iI] - 1].x));
						break;
					case TSIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].x + pa[nvtx[1][sosedb[i].iI] - 1].x));
						break;
					case BSIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].x + pa[nvtx[1][sosedb[i].iI] - 1].x));
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// запись y
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (sosedb[i].Norm) {// определим внутреннюю нормаль к границе
					case ESIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].y + pa[nvtx[2][sosedb[i].iI] - 1].y));
						break;
					case WSIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].y + pa[nvtx[2][sosedb[i].iI] - 1].y));
						break;
					case NSIDE : fprintf(fp, "%+.16f ", pa[nvtx[0][sosedb[i].iI] - 1].y);
						break;
					case SSIDE: fprintf(fp, "%+.16f ", pa[nvtx[2][sosedb[i].iI] - 1].y);
						break;
					case TSIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].y + pa[nvtx[2][sosedb[i].iI] - 1].y));
						break;
					case BSIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].y + pa[nvtx[2][sosedb[i].iI] - 1].y));
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// запись z
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (sosedb[i].Norm) {// определим внутреннюю нормаль к границе
					case ESIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].z + pa[nvtx[4][sosedb[i].iI] - 1].z));
						break;
					case WSIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].z + pa[nvtx[4][sosedb[i].iI] - 1].z));
						break;
					case NSIDE : fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].z + pa[nvtx[4][sosedb[i].iI] - 1].z));
						break;
					case SSIDE:fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][sosedb[i].iI] - 1].z + pa[nvtx[4][sosedb[i].iI] - 1].z));
						break;
					case TSIDE: fprintf(fp, "%+.16f ", pa[nvtx[0][sosedb[i].iI] - 1].z);
						break;
					case BSIDE: fprintf(fp, "%+.16f ", pa[nvtx[4][sosedb[i].iI] - 1].z);
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			fclose(fp); // закрытие файла
			printf("export tecplot temperature part1 is successfully written...OK.\n");

		}

#ifdef MINGW_COMPILLER
		err = 0;
		fp=fopen64("ALICEFLOW0_06_temp_part3.txt", "w");
		if (fp == nullptr) err = 1;
#else
		err = fopen_s(&fp, "ALICEFLOW0_06_temp_part3.txt", "w");
#endif
		if ((err) != 0) {
			printf("Create File temp part3 Error\n");
			//getchar();
			system("pause");
		}
		else {

			if (fp != nullptr) {

				integer i = 0;
				// запись информации о разностной сетке
				for (i = 0; i < ncell; i++) {
#if doubleintprecision == 1
					fprintf(fp, "%lld %lld %lld %lld ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#else
					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#endif
					
				}

				fclose(fp); // закрытие файла
				printf("export tecplot temperature part3 is successfully written...OK.\n");
			}
		}
	}

} // exporttecplotxy360T_3D_part1and3

// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2binary(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint)
{
    // ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage=false;

	FILE *fp=NULL;
    FILE *fp1=NULL; // часть 1 или 3
	errno_t err=0;
#ifdef MINGW_COMPILLER
	fp=fopen64("ALICEFLOW0_07_temp.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");
#endif
	
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	
	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
	if ((err) != 0) {
		printf("Create File temp Error in exporttecplotxy360T_3D_part2binary in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {
        
        int c; // читаемый символ
		integer ivarexport=1; // по умолчанию только поле температур:
		integer i=0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1=fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1 == NULL) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		char symbol;
		doublereal fnumber;
		if (bOk)
		{

        
			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i=0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport=3; // считаем что температура всегда активна
				}
			}
            
			if (ivarexport==1) {
				// запись заголовка
		        //fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");
				char title[] = "TITLE = ";
				for (i = 0; i <= 7; i++) {
					symbol = title[i];
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				 symbol = '\"';
				fwrite(&symbol, sizeof(char), 1, fp);
				char tit7[] = "ALICEFLOW0_07";
				for (i = 0; i <= 12; i++) {
					symbol = tit7[i];
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				 symbol = '\"';
				fwrite(&symbol, sizeof(char), 1, fp);
				 symbol = '\n';
				fwrite(&symbol, sizeof(char), 1, fp);


		         // запись имён переменных
		        // fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
				 char variablestit[] = "VARIABLES = x, y, z, Temp, Lam";
				 //31
				 for (i = 0; i <= 29; i++) {
					 symbol = variablestit[i];
					 fwrite(&symbol, sizeof(char), 1, fp);
				 }
				  symbol = '\n';
				 fwrite(&symbol, sizeof(char), 1, fp);


				 // запись информации о зонах
				 if (bextendedprint) {
#if doubleintprecision == 1
					 //fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm+t.maxbound, ncell);
#else
					 //fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm+t.maxbound, ncell);
#endif

					
					 char tit3[] = "ZONE T=";//20
					 for (i = 0; i <= 6; i++) {
						 symbol = tit3[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					  symbol = '\"';
					 fwrite(&symbol, sizeof(char), 1, fp);
					 char tit12[] = "Rampant";
					 for (i = 0; i <= 6; i++) {
						 symbol = tit12[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					  symbol = '\"';
					 fwrite(&symbol, sizeof(char), 1, fp);
					 char tit13[] = ", N=";
					 for (i = 0; i <= 3; i++) {
						 symbol = tit13[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					 integer number = maxelm + t.maxbound;
					 fwrite(&number, sizeof(integer), 1, fp);
					 char tit4[] = ", E=";//4
					 for (i = 0; i <= 3; i++) {
						 symbol = tit4[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					 fwrite(&ncell, sizeof(integer), 1, fp);
					 char tit5[] = ", ET=BRICK, F=FEBLOCK";//23
					 for (i = 0; i <= 20; i++) {
						 symbol = tit5[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					  symbol = '\n';
					 fwrite(&symbol, sizeof(char), 1, fp);
					 // symbol = '\n';
					 //fwrite(&symbol, sizeof(char), 1, fp);
				 }
				 else {
#if doubleintprecision == 1
					 // fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
					 // fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif

                  
					 char tit3[] = "ZONE T=";//20
					 for (i = 0; i <= 6; i++) {
						 symbol = tit3[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					  symbol = '\"';
					 fwrite(&symbol, sizeof(char), 1, fp);
					 char tit12[] = "Rampant";
					 for (i = 0; i <= 6; i++) {
						 symbol = tit12[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					  symbol = '\"';
					 fwrite(&symbol, sizeof(char), 1, fp);
					 char tit13[] = ", N=";
					 for (i = 0; i <= 3; i++) {
						 symbol = tit13[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					 integer number = maxelm;
					 fwrite(&number, sizeof(integer), 1, fp);
					 char tit4[] = ", E=";//4
					 for (i = 0; i <= 3; i++) {
						 symbol = tit4[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					 fwrite(&ncell, sizeof(integer), 1, fp);
					 char tit5[] = ", ET=BRICK, F=FEBLOCK";//23
					 for (i = 0; i <= 20; i++) {
						 symbol = tit5[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					  symbol = '\n';
					 fwrite(&symbol, sizeof(char), 1, fp);
					//  symbol = '\n';
					 //fwrite(&symbol, sizeof(char), 1, fp);

				 }
		
				 if (bvery_big_memory) {
					 // extended printeger не предусмотрено.

					 // запись x
					 for (i = 0; i < t.database.maxelm; i++) {
						 //fprintf(fp, "%+.16f ", t.database.x[i]);
						 fnumber = t.database.x[i];
						 fwrite(&fnumber, sizeof(doublereal), 1, fp);
						 //symbol = ' ';
						 //fwrite(&symbol, sizeof(char), 1, fp);
						 //if (i % 10 == 0) {
							 //fprintf(fp, "\n");
							//char symbol = '\n';
							 //fwrite(&symbol, sizeof(char), 1, fp);
						 //}
					 }
					 // запись y
					 for (i = 0; i < t.database.maxelm; i++) {
						 //fprintf(fp, "%+.16f ", t.database.y[i]);
						 fnumber = t.database.y[i];
						 fwrite(&fnumber, sizeof(doublereal), 1, fp);
						 //symbol = ' ';
						 //fwrite(&symbol, sizeof(char), 1, fp);
						 //if (i % 10 == 0) {
							 //fprintf(fp, "\n");
							//char symbol = '\n';
							 //fwrite(&symbol, sizeof(char), 1, fp);
						 //}
					 }
					 // запись z
					 for (i = 0; i < t.database.maxelm; i++) {
						 //fprintf(fp, "%+.16f ", t.database.z[i]);
						 fnumber = t.database.z[i];
						 fwrite(&fnumber, sizeof(doublereal), 1, fp);
						 //symbol = ' ';
						 //fwrite(&symbol, sizeof(char), 1, fp);
						 //if (i % 10 == 0) {
							 //fprintf(fp, "\n");
							// char symbol = '\n';
							 //fwrite(&symbol, sizeof(char), 1, fp);
						// }
					 }
				 }
				 else {
					 while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				 }
			}
			else if (ivarexport==3) {
				// запись заголовка
		        fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// Полный набор искомых величин и теплопередача и гидродинамика:
				fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z\n");
					
#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#endif

			
				if (bvery_big_memory) {
					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						//fprintf(fp, "%+.16f ", t.database.x[i]);
						fnumber = t.database.x[i];
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
						if (i % 10 == 0) {
							//fprintf(fp, "\n");
							char symbol = '\n';
							fwrite(&symbol, sizeof(char), 1, fp);
						}
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						//fprintf(fp, "%+.16f ", t.database.y[i]);
						fnumber = t.database.y[i];
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
						if (i % 10 == 0) {
							//fprintf(fp, "\n");
							char symbol = '\n';
							fwrite(&symbol, sizeof(char), 1, fp);
						}
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						//fprintf(fp, "%+.16f ", t.database.z[i]);
						fnumber = t.database.z[i];
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
						if (i % 10 == 0) {
							//fprintf(fp, "\n");
							char symbol = '\n';
							fwrite(&symbol, sizeof(char), 1, fp);
						}
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части
        
        // Запись поля температур производится всегда.
		

		// запись температуры
		for (i=0;i<maxelm; i++) {
			//fprintf(fp, "%+.16f ", t.potent[i]);
			fnumber = t.potent[i];
			fwrite(&fnumber, sizeof(doublereal), 1, fp);
			//symbol = ' ';
			//fwrite(&symbol, sizeof(char), 1, fp);
			//if (i % 10 == 0) {
				//fprintf(fp, "\n");
			//	symbol = '\n';
			//	fwrite(&symbol, sizeof(char), 1, fp);
			//}
		}

		if (bextendedprint) {
			for (i=maxelm;i<maxelm+t.maxbound; i++) {
			  //  fprintf(fp, "%+.16f ", t.potent[i]);
				fnumber = t.potent[i];
				fwrite(&fnumber, sizeof(doublereal), 1, fp);
				//symbol = ' ';
				//fwrite(&symbol, sizeof(char), 1, fp);
				//if ((i + maxelm) % 10 == 0) {
					//fprintf(fp, "\n");
					//symbol = '\n';
					//fwrite(&symbol, sizeof(char), 1, fp);
				//}
			}
		}
		

		//fprintf(fp, "\n");
		//symbol = '\n';
		//fwrite(&symbol, sizeof(char), 1, fp);

		// Lam
		for (i=0;i<maxelm; i++) {
			//fprintf(fp, "%+.16f ", t.prop[LAM][i]);
			fnumber = t.prop[LAM][i];
			fwrite(&fnumber, sizeof(doublereal), 1, fp);
			//symbol = ' ';
			//fwrite(&symbol, sizeof(char), 1, fp);
			//if (i % 10 == 0) {
				//fprintf(fp, "\n");
			//	symbol = '\n';
			//	fwrite(&symbol, sizeof(char), 1, fp);
			//}
		}

		if (bextendedprint) {
			for (i=0;i<t.maxbound; i++) {
			    // fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
				 fnumber = t.prop_b[LAM][i];
				 fwrite(&fnumber, sizeof(doublereal), 1, fp);
				 //symbol = ' ';
				 //fwrite(&symbol, sizeof(char), 1, fp);
				 //if ((i + maxelm) % 10 == 0) {
					 //fprintf(fp, "\n");
					// symbol = '\n';
					 //fwrite(&symbol, sizeof(char), 1, fp);
				// }
		    }
		}

		//fprintf(fp, "\n");
		//symbol = '\n';
		//fwrite(&symbol, sizeof(char), 1, fp);

		// Запись гидродинамических величин если необходимо:
		if (ivarexport==3) {
			char symbol;
			doublereal fnumber;
			// Speed
            for (i=0;i<maxelm; i++) {
				if (t.ptr[1][i]>-1) {
					doublereal svx=f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy=f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz=f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					//fprintf(fp, "%+.16f ",sqrt(svx+svy+svz)); 
					fnumber = sqrt(svx + svy + svz);
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			
			if (bextendedprint) {
				// Speed
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
				    	doublereal svx=f[idfluid].potent[VX][i+maxelm]*f[idfluid].potent[VX][i+maxelm];
					    doublereal svy=f[idfluid].potent[VY][i+maxelm]*f[idfluid].potent[VY][i+maxelm];
					    doublereal svz=f[idfluid].potent[VZ][i+maxelm]*f[idfluid].potent[VZ][i+maxelm];
					    //fprintf(fp, "%+.16f ",sqrt(svx+svy+svz));
						fnumber = sqrt(svx + svy + svz);
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
						if ((i + maxelm) % 10 == 0) {
							//fprintf(fp, "\n");
							symbol = '\n';
							fwrite(&symbol, sizeof(char), 1, fp);
						}
		        }
			}
			

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// Pressure
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
				   fnumber = f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// Pressure
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       //fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i+maxelm]); // PRESSURE
					   fnumber = f[idfluid].potent[PRESS][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// PAM
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   //fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				   fnumber = f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// PAM
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i+maxelm]); // PRESSURE
					   fnumber = f[idfluid].potent[PAM][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// VX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                 //  fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// VX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i+maxelm]); // VX
					   fnumber = f[idfluid].potent[VX][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// VY
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// VY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i+maxelm]); // VY
					   fnumber = f[idfluid].potent[VY][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);


			// VZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// VZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i+maxelm]); // VZ
					   fnumber = f[idfluid].potent[VZ][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);
			
            // Rho
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					fnumber = f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]];
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", t.prop[RHO][i]);
					fnumber = t.prop[RHO][i];
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// Rho
				for (i=0;i<f[0].maxbound; i++) {
					//fprintf(fp, "%+.16f ",f[0].prop_b[RHO][i]);
					fnumber = f[0].prop_b[RHO][i];
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
					if ((i + maxelm) % 10 == 0) {
						//fprintf(fp, "\n");
						symbol = '\n';
						fwrite(&symbol, sizeof(char), 1, fp);
					}
				}
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// Mu
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//--->//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
					fnumber = f[t.ptr[1][i]].prop[MU][t.ptr[0][i]];
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// Mu
				for (i=0;i<f[0].maxbound; i++) {
					//fprintf(fp, "%+.16f ",f[0].prop_b[MU][i]);
					fnumber = f[0].prop_b[MU][i];
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
					if ((i + maxelm) % 10 == 0) {
						//fprintf(fp, "\n");
						symbol = '\n';
						fwrite(&symbol, sizeof(char), 1, fp);
					}
				}
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// Mut // Турбулентная динамическая вязкость
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                 //  fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// MUT
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i+maxelm]); // MUT
					   fnumber = f[idfluid].potent[MUT][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   //fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime==ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime==SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES)||
						(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
						fnumber = f[t.ptr[1][i]].rdistWall[t.ptr[0][i]];
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
					}
					else {
						//fprintf(fp, "%+.16f ", 0.0);
						fnumber = 0.0;
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
					}
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime==ZEROEQMOD) || 
						(f[0].iflowregime==SMAGORINSKY)||
						(f[0].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == RANS_MENTER_SST)||
						(f[0].iflowregime == RANS_STANDART_K_EPS)) {
                      // fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i+maxelm]); // Distance_Wall
					   fnumber = f[idfluid].rdistWall[i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					}
					else {
						//fprintf(fp, "%+.16f ", 0.0);
						fnumber = 0.0;
						fwrite(&fnumber, sizeof(doublereal), 1, fp);
						symbol = ' ';
						fwrite(&symbol, sizeof(char), 1, fp);
					}
					if ((i + maxelm) % 10 == 0) {
						//fprintf(fp, "\n");
						symbol = '\n';
						fwrite(&symbol, sizeof(char), 1, fp);
					}
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			

			// Curl // Завихрённость - модуль ротора скорости.
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                 //  fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				   fnumber = f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// Curl
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i+maxelm]); // Curl
					   fnumber = f[idfluid].potent[CURL][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADXVX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       //fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i+maxelm]); // GRADXVX
					   fnumber = f[idfluid].potent[GRADXVX][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADYVX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   //fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADYVX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                     //  fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i+maxelm]); // GRADYVX
					   fnumber = f[idfluid].potent[GRADYVX][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADZVX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                 //  fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADZVX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i+maxelm]); // GRADZVX
					   fnumber = f[idfluid].potent[GRADZVX][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADXVY
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]); 
				   fnumber = f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADXVY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i+maxelm]); // GRADXVY
					   fnumber = f[idfluid].potent[GRADXVY][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADYVY
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
				    //	fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADYVY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i+maxelm]); // GRADYVY
					   fnumber = f[idfluid].potent[GRADYVY][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADZVY
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                 //  fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADZVY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       //fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i+maxelm]); // GRADZVY
					   fnumber = f[idfluid].potent[GRADZVY][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADXVZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i+maxelm]); // GRADXVZ
					   fnumber = f[idfluid].potent[GRADXVZ][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADYVZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   //fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       //fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i+maxelm]); // GRADYVZ
					   fnumber = f[idfluid].potent[GRADYVZ][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						   //fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// GRADZVZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				}
				else {
					//fprintf(fp, "%+.16f ", 0.0);
					fnumber = 0.0;
					fwrite(&fnumber, sizeof(doublereal), 1, fp);
					symbol = ' ';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       //fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i+maxelm]); // GRADZVZ
					   fnumber = f[idfluid].potent[GRADZVZ][i + maxelm];
					   fwrite(&fnumber, sizeof(doublereal), 1, fp);
					   symbol = ' ';
					   fwrite(&symbol, sizeof(char), 1, fp);
					   if ((i + maxelm) % 10 == 0) {
						  // fprintf(fp, "\n");
						   symbol = '\n';
						   fwrite(&symbol, sizeof(char), 1, fp);
					   }
		        }
			}

		   // fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			doublereal *Tx=nullptr;
            doublereal *Ty=nullptr;
            doublereal *Tz=nullptr;
			Tx=new doublereal[t.maxelm+t.maxbound];
			Ty=new doublereal[t.maxelm+t.maxbound];
			Tz=new doublereal[t.maxelm+t.maxbound];

			// инициализация нулём.
			for (i=0; i<t.maxelm+t.maxbound; i++) {
				Tx[i]=0.0;
				Ty[i]=0.0;
				Tz[i]=0.0;
			}

			// нахождение градиентов.
			for (i=0; i<t.maxelm; i++) {
				// Только внутренние узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, false, 
					t.sosedb,  Tx, Ty, Tz, t.ilevel_alice);
			}

			for (i=0; i<t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, true, 
					t.sosedb,Tx, Ty, Tz, t.ilevel_alice);
			}
			

			// Сохранение в файл.

			// Heat Flux X
		    for (i=0;i<maxelm; i++) {
			    //fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Tx[i]);
				fnumber = -t.prop[LAM][i] * Tx[i];
				fwrite(&fnumber, sizeof(doublereal), 1, fp);
				symbol = ' ';
				fwrite(&symbol, sizeof(char), 1, fp);
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			      // fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Tx[i+maxelm]);
				   fnumber = -t.prop_b[LAM][i] * Tx[i + maxelm];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				   if ((i + maxelm) % 10 == 0) {
					   //fprintf(fp, "\n");
					   symbol = '\n';
					   fwrite(&symbol, sizeof(char), 1, fp);
				   }
		       }
		    }

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// Heat Flux Y
		    for (i=0;i<maxelm; i++) {
			    //fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Ty[i]);
				fnumber = -t.prop[LAM][i] * Ty[i];
				fwrite(&fnumber, sizeof(doublereal), 1, fp);
				symbol = ' ';
				fwrite(&symbol, sizeof(char), 1, fp);
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       //fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Ty[i+maxelm]);
				   fnumber = -t.prop_b[LAM][i] * Ty[i + maxelm];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				   if ((i + maxelm) % 10 == 0) {
					  // fprintf(fp, "\n");
					   symbol = '\n';
					   fwrite(&symbol, sizeof(char), 1, fp);
				   }
		       }
		    }

		    //fprintf(fp, "\n");
			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);

			// Heat Flux Z
		    for (i=0;i<maxelm; i++) {
			    //fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Tz[i]);
				fnumber = -t.prop[LAM][i] * Tz[i];
				fwrite(&fnumber, sizeof(doublereal), 1, fp);
				symbol = ' ';
				fwrite(&symbol, sizeof(char), 1, fp);
				if (i % 10 == 0) {
					//fprintf(fp, "\n");
					symbol = '\n';
					fwrite(&symbol, sizeof(char), 1, fp);
				}
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       //fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Tz[i+maxelm]);
				   fnumber = -t.prop_b[LAM][i] * Tz[i + maxelm];
				   fwrite(&fnumber, sizeof(doublereal), 1, fp);
				   symbol = ' ';
				   fwrite(&symbol, sizeof(char), 1, fp);
				   if ((i + maxelm) % 10 == 0) {
					  // fprintf(fp, "\n");
					  symbol = '\n';
					  fwrite(&symbol, sizeof(char), 1, fp);
				   }
		       }
		    }

			symbol = '\n';
			fwrite(&symbol, sizeof(char), 1, fp);
		   // fprintf(fp, "\n");

			// Освобождение оперативной памяти.
			if (Tx != nullptr) {
				delete[] Tx;
			}
			if (Ty != nullptr) {
				delete[] Ty;
			}
			if (Tz != nullptr) {
				delete[] Tz;
			}

		}

		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
#if doubleintprecision == 1
				//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
				//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);

#else
				//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
				//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);

#endif
				integer number = t.database.nvtxcell[0][i];
				fwrite(&number, sizeof(integer), 1, fp);
				char symbol = ' ';
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[1][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[2][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[3][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[4][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[5][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[6][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = t.database.nvtxcell[7][i];
				fwrite(&number, sizeof(integer), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				symbol = '\n';
				fwrite(&symbol, sizeof(char), 1, fp);
			}
		}
		/*
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}
		*/

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}

// 30.09.2019 решено сохранять также и турбулентную вязкость в случае турбулентного течения.
void save_velocity_for_init(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count) {


	errno_t err_inicialization_data=0;
	FILE* fp_inicialization_data=NULL;
#ifdef MINGW_COMPILLER
	fp_inicialization_data=fopen64("load.txt", "w");
	if (fp_inicialization_data==NULL) err_inicialization_data = 1;
#else
	err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "w");
#endif

	
	if (err_inicialization_data != 0) {
		// открытие неудачно или файл отсутствует.
		printf("Create File load.txt Error\n");
		//getchar();
		system("pause");
	}
	else {

		if (b_on_adaptive_local_refinement_mesh) {
			// Гидродинамика была посчитана на АЛИС сетке.


			// Работает только для одной FLUID зоны.
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
			fprintf(fp_inicialization_data, "%lld\n %lld\n %lld\n", f[0].maxnod, f[0].maxelm, f[0].iflowregime);
#else
			fprintf(fp_inicialization_data, "%d\n %d\n %d\n", f[0].maxnod, f[0].maxelm, f[0].iflowregime);
#endif
			// X coordinate
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", f[0].pa[i].x);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");

			// Y coordinate
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", f[0].pa[i].y);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");

			// Z coordinate
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", f[0].pa[i].z);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");
#else
#if doubleintprecision == 1
			fprintf_s(fp_inicialization_data, "%lld\n %lld\n %lld\n", f[0].maxnod, f[0].maxelm, f[0].iflowregime);
#else
			fprintf_s(fp_inicialization_data, "%d\n %d\n %d\n", f[0].maxnod, f[0].maxelm, f[0].iflowregime);
#endif
			// X coordinate
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", f[0].pa[i].x);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
			fprintf_s(fp_inicialization_data, "\n");

			// Y coordinate
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", f[0].pa[i].y);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
			}
			fprintf_s(fp_inicialization_data, "\n");

			// Z coordinate
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", f[0].pa[i].z);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
			}
			fprintf_s(fp_inicialization_data, "\n");
#endif



			// Переносим компоненты скорости из центров ячеек в узлы сетки.
			// Простая и надежная интерполяция.
			doublereal* Ux = nullptr;
			if (f[0].maxnod > 0) {
				Ux = new doublereal[f[0].maxnod];
			}
			doublereal* Uy = nullptr;
			if (f[0].maxnod > 0) {
				Uy = new doublereal[f[0].maxnod];
			}
			doublereal* Uz = nullptr;
			if (f[0].maxnod > 0) {
				Uz = new doublereal[f[0].maxnod];
			}
			doublereal* mut = nullptr;
			if (f[0].maxnod > 0) {
				mut = new doublereal[f[0].maxnod];
			}
			doublereal* vol = nullptr;
			if (f[0].maxnod > 0) {
				vol = new doublereal[f[0].maxnod];
			}
			//doublereal* vesaX = new doublereal[f[0].maxnod];
			//doublereal* vesaY = new doublereal[f[0].maxnod];
			//doublereal* vesaZ = new doublereal[f[0].maxnod];

			for (integer i = 0; i < f[0].maxnod; i++) {
				Ux[i] = 0.0;
				Uy[i] = 0.0;
				Uz[i] = 0.0;
				mut[i] = 0.0;
				vol[i] = 0.0;
				//vesaX[i] = 0.0;
				//vesaY[i] = 0.0;
				//vesaZ[i] = 0.0;
			}

			doublereal SpeedMax = -1.0e60;
			doublereal SpeedMin = 1.0e60;

			for (integer i = 0; i < f[0].maxelm; i++) {
				integer inode1 = 0, inode2 = 0, inode3 = 0, inode4 = 0, inode5 = 0, inode6 = 0, inode7 = 0, inode8 = 0;
				inode1 = f[0].nvtx[0][i] - 1;
				inode2 = f[0].nvtx[1][i] - 1;
				inode3 = f[0].nvtx[2][i] - 1;
				inode4 = f[0].nvtx[3][i] - 1;
				inode5 = f[0].nvtx[4][i] - 1;
				inode6 = f[0].nvtx[5][i] - 1;
				inode7 = f[0].nvtx[6][i] - 1;
				inode8 = f[0].nvtx[7][i] - 1;
				if ((inode1 < 0) || (inode2 < 0) || (inode3 < 0) || (inode4 < 0) || (inode5 < 0) || (inode6 < 0) || (inode7 < 0) || (inode8 < 0)) {
					printf("ERROR: NEGATIVE NVTX in function save_velocity_for_init in module my_export_tecplot3.c\n");
					system("PAUSE");
					exit(1);
				}

				doublereal Speed = sqrt((f[0].potent[VX][i])*(f[0].potent[VX][i])+(f[0].potent[VY][i])*(f[0].potent[VY][i])+(f[0].potent[VZ][i])*(f[0].potent[VZ][i]));
				if (Speed > SpeedMax) SpeedMax = Speed;
				if (Speed < SpeedMin) SpeedMin = Speed;

				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, f[0].nvtx, f[0].pa, dx, dy, dz);
				
				if (Ux != nullptr) {
					Ux[inode1] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode2] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode3] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode4] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode5] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode6] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode7] += 0.125*dx*dy*dz*f[0].potent[VX][i];
					Ux[inode8] += 0.125*dx*dy*dz*f[0].potent[VX][i];
				}

				if (Uy != nullptr) {
					Uy[inode1] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode2] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode3] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode4] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode5] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode6] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode7] += 0.125*dx*dy*dz*f[0].potent[VY][i];
					Uy[inode8] += 0.125*dx*dy*dz*f[0].potent[VY][i];
				}

				if (Uz != nullptr) {
					Uz[inode1] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode2] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode3] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode4] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode5] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode6] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode7] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
					Uz[inode8] += 0.125*dx*dy*dz*f[0].potent[VZ][i];
				}

				if (mut != nullptr) {
					mut[inode1] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode2] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode3] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode4] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode5] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode6] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode7] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
					mut[inode8] += 0.125*dx*dy*dz*f[0].potent[MUT][i];
				}
				/*
				Ux[inode1] += (1.0/(0.5*dx))*f[0].potent[VX][i];
				Ux[inode2] += (1.0/(0.5*dx))*f[0].potent[VX][i];
				Ux[inode3] += (1.0/(0.5*dx))*f[0].potent[VX][i];
				Ux[inode4] += (1.0 / (0.5*dx))*f[0].potent[VX][i];
				Ux[inode5] += (1.0 / (0.5*dx))*f[0].potent[VX][i];
				Ux[inode6] += (1.0 / (0.5*dx))*f[0].potent[VX][i];
				Ux[inode7] += (1.0 / (0.5*dx))*f[0].potent[VX][i];
				Ux[inode8] += (1.0/(0.5*dx))*f[0].potent[VX][i];

				vesaX[inode1] += (1.0 / (0.5*dx));
				vesaX[inode2] += (1.0 / (0.5*dx));
				vesaX[inode3] += (1.0 / (0.5*dx));
				vesaX[inode4] += (1.0 / (0.5*dx));
				vesaX[inode5] += (1.0 / (0.5*dx));
				vesaX[inode6] += (1.0 / (0.5*dx));
				vesaX[inode7] += (1.0 / (0.5*dx));
				vesaX[inode8] += (1.0 / (0.5*dx));

				Uy[inode1] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode2] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode3] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode4] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode5] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode6] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode7] += (1.0 / (0.5*dy))*f[0].potent[VY][i];
				Uy[inode8] += (1.0 / (0.5*dy))*f[0].potent[VY][i];

				vesaY[inode1] += (1.0 / (0.5*dy));
				vesaY[inode2] += (1.0 / (0.5*dy));
				vesaY[inode3] += (1.0 / (0.5*dy));
				vesaY[inode4] += (1.0 / (0.5*dy));
				vesaY[inode5] += (1.0 / (0.5*dy));
				vesaY[inode6] += (1.0 / (0.5*dy));
				vesaY[inode7] += (1.0 / (0.5*dy));
				vesaY[inode8] += (1.0 / (0.5*dy));


				Uz[inode1] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode2] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode3] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode4] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode5] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode6] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode7] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];
				Uz[inode8] += (1.0 / (0.5*dz))*f[0].potent[VZ][i];

				vesaZ[inode1] += (1.0 / (0.5*dz));
				vesaZ[inode2] += (1.0 / (0.5*dz));
				vesaZ[inode3] += (1.0 / (0.5*dz));
				vesaZ[inode4] += (1.0 / (0.5*dz));
				vesaZ[inode5] += (1.0 / (0.5*dz));
				vesaZ[inode6] += (1.0 / (0.5*dz));
				vesaZ[inode7] += (1.0 / (0.5*dz));
				vesaZ[inode8] += (1.0 / (0.5*dz));
				*/
				if (vol != nullptr) {
					vol[inode1] += 0.125*dx*dy*dz;
					vol[inode2] += 0.125*dx*dy*dz;
					vol[inode3] += 0.125*dx*dy*dz;
					vol[inode4] += 0.125*dx*dy*dz;
					vol[inode5] += 0.125*dx*dy*dz;
					vol[inode6] += 0.125*dx*dy*dz;
					vol[inode7] += 0.125*dx*dy*dz;
					vol[inode8] += 0.125*dx*dy*dz;
				}
			}

			printf("Velocity statistics:\n");
			printf("Calculate vel: SpeedMax = %e SpeedMin = %e\n", SpeedMax, SpeedMin);
			doublereal Speed_coef = SpeedMax;
			SpeedMax = -1.0e60;
			SpeedMin = 1.0e60;

			// Средневзвешенное по объёму значение компоненты скорости.
			for (integer i = 0; i < f[0].maxnod; i++) {
				if (fabs(vol[i]) < 1.0e-40) {
					Ux[i] = 0.0;
					Uy[i] = 0.0;
					Uz[i] = 0.0;
					mut[i] = 0.0;
				}
				else {
			        Ux[i] /= vol[i];
					Uy[i] /= vol[i];
					Uz[i] /= vol[i];
					mut[i] /= vol[i];
					//Ux[i] /= vesaX[i];
					//Uy[i] /= vesaY[i];
					//Uz[i] /= vesaZ[i];
					doublereal Speed = sqrt((Ux[i])*(Ux[i])+(Uy[i])*(Uy[i])+(Uz[i])*(Uz[i]));
					if (Speed > SpeedMax) SpeedMax = Speed;
					if (Speed < SpeedMin) SpeedMin = Speed;
				}
			}

			//Домножение на корректирующий коэффициент приближает поле температур ближе 
			// к истинному. 21.03.2019
			Speed_coef /= SpeedMax;// Коэффициент больше 1. (например 1.6).
			for (integer i = 0; i < f[0].maxnod; i++) {
				Ux[i] *= Speed_coef;
				Uy[i] *= Speed_coef;
				Uz[i] *= Speed_coef;
			}


			printf("Interpolation vel: SpeedMax = %e SpeedMin = %e\n", SpeedMax, SpeedMin);

			delete[] vol;
			vol = nullptr;

			//delete[] vesaX;
			//vesaX = nullptr;
			//delete[] vesaY;
			//vesaY = nullptr;
			//delete[] vesaZ;
			//vesaZ = nullptr;

#ifdef MINGW_COMPILLER
            // UX speed component
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", Ux[i]);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");

			// UY speed component
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", Uy[i]);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");

			// UZ speed component
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", Uz[i]);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");

			// турбулентная динамическая вязкость.
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf(fp_inicialization_data, "%e ", mut[i]);
				if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
			}
			fprintf(fp_inicialization_data, "\n");
#else
            // UX speed component
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", Ux[i]);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
			}
			fprintf_s(fp_inicialization_data, "\n");

			// UY speed component
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", Uy[i]);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
			}
			fprintf_s(fp_inicialization_data, "\n");

			// UZ speed component
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", Uz[i]);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
			}
			fprintf_s(fp_inicialization_data, "\n");

			// турбулентная динамическая вязкость.
			for (integer i = 0; i < f[0].maxnod; i++) {
				fprintf_s(fp_inicialization_data, "%e ", mut[i]);
				if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
			}
			fprintf_s(fp_inicialization_data, "\n");
#endif
			

			delete[] Ux;
			Ux = nullptr;
			delete[] Uy;
			Uy = nullptr;
			delete[] Uz;
			Uz = nullptr;
			delete[] mut;
			mut = nullptr;

			for (integer i = 0; i < f[0].maxelm; i++) {
				integer inode1=0, inode2=0, inode3=0, inode4=0, inode5=0, inode6=0, inode7=0, inode8=0;
				inode1 = f[0].nvtx[0][i] - 1;
				inode2 = f[0].nvtx[1][i] - 1;
				inode3 = f[0].nvtx[2][i] - 1;
				inode4 = f[0].nvtx[3][i] - 1;
				inode5 = f[0].nvtx[4][i] - 1;
				inode6 = f[0].nvtx[5][i] - 1;
				inode7 = f[0].nvtx[6][i] - 1;
				inode8 = f[0].nvtx[7][i] - 1;
				if ((inode1 < 0) || (inode2 < 0) || (inode3 < 0) || (inode4 < 0) || (inode5 < 0) || (inode6 < 0) || (inode7 < 0) || (inode8 < 0)) {
					printf("ERROR: NEGATIVE NVTX in function save_velocity_for_init in module my_export_tecplot3.c\n");
					system("PAUSE");
					exit(1);
				}
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
				fprintf(fp_inicialization_data, "%lld %lld %lld %lld %lld %lld %lld %lld\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#else
				fprintf(fp_inicialization_data, "%d %d %d %d %d %d %d %d\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#endif
#else
#if doubleintprecision == 1
				fprintf_s(fp_inicialization_data, "%lld %lld %lld %lld %lld %lld %lld %lld\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#else
				fprintf_s(fp_inicialization_data, "%d %d %d %d %d %d %d %d\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#endif
#endif

			}

			// Экспорт на АЛИС сетке завершен.

		}
		else {

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
		fprintf(fp_inicialization_data, "%lld\n %lld\n %lld\n", maxelm, ncell, f[0].iflowregime);
#else
		fprintf(fp_inicialization_data, "%d\n %d\n %d\n", maxelm, ncell, f[0].iflowregime);
#endif
		// X
		for (integer i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf(fp_inicialization_data, "%e ", p.x);
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");
		// Y
		for (integer i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf(fp_inicialization_data, "%e ", p.y);
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");
		// Z
		for (integer i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf(fp_inicialization_data, "%e ", p.z);
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");

		// VX
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf(fp_inicialization_data, "%e ", f[0].potent[VX][t.ptr[0][i]]);
			}
			else {
				fprintf(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");

		// VY
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf(fp_inicialization_data, "%e ", f[0].potent[VY][t.ptr[0][i]]);
			}
			else {
				fprintf(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");

		// VZ
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf(fp_inicialization_data, "%e ", f[0].potent[VZ][t.ptr[0][i]]);
			}
			else {
				fprintf(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");

		// MUT
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf(fp_inicialization_data, "%e ", f[0].potent[MUT][t.ptr[0][i]]);
			}
			else {
				fprintf(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf(fp_inicialization_data, "\n");
		}
		fprintf(fp_inicialization_data, "\n");
#else
#if doubleintprecision == 1
		fprintf_s(fp_inicialization_data, "%lld\n %lld\n %lld\n", maxelm, ncell, f[0].iflowregime);
#else
		fprintf_s(fp_inicialization_data, "%d\n %d\n %d\n", maxelm, ncell, f[0].iflowregime);
#endif

		// X
		for (integer i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf_s(fp_inicialization_data, "%e ", p.x);
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		// Y
		for (integer i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf_s(fp_inicialization_data, "%e ", p.y);
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		// Z
		for (integer i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf_s(fp_inicialization_data, "%e ", p.z);
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		
		// VX
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr != nullptr) {
				if (t.ptr[1][i] > -1) {
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VX][t.ptr[0][i]]);
				}
				else {
					fprintf_s(fp_inicialization_data, "%e ", 0.0);
				}
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		// VY
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr != nullptr) {
				if (t.ptr[1][i] > -1) {
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VY][t.ptr[0][i]]);
				}
				else {
					fprintf_s(fp_inicialization_data, "%e ", 0.0);
				}
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		// VZ
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr != nullptr) {
				if (t.ptr[1][i] > -1) {
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VZ][t.ptr[0][i]]);
				}
				else {
					fprintf_s(fp_inicialization_data, "%e ", 0.0);
				}
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		// MUT
		for (integer i = 0; i < maxelm; i++) {
			if (t.ptr != nullptr) {
				if (t.ptr[1][i] > -1) {
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[MUT][t.ptr[0][i]]);
				}
				else {
					fprintf_s(fp_inicialization_data, "%e ", 0.0);
				}
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");

#endif			

			for (integer i = 0; i < ncell; i++) {
				integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = t.database.nvtxcell[0][i] - 1;
				inode2 = t.database.nvtxcell[1][i] - 1;
				inode3 = t.database.nvtxcell[2][i] - 1;
				inode4 = t.database.nvtxcell[3][i] - 1;
				inode5 = t.database.nvtxcell[4][i] - 1;
				inode6 = t.database.nvtxcell[5][i] - 1;
				inode7 = t.database.nvtxcell[6][i] - 1;
				inode8 = t.database.nvtxcell[7][i] - 1;
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
				fprintf(fp_inicialization_data, "%lld %lld %lld %lld %lld %lld %lld %lld\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#else
				fprintf(fp_inicialization_data, "%d %d %d %d %d %d %d %d\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#endif
#else
#if doubleintprecision == 1
				fprintf_s(fp_inicialization_data, "%lld %lld %lld %lld %lld %lld %lld %lld\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#else
				fprintf_s(fp_inicialization_data, "%d %d %d %d %d %d %d %d\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#endif
#endif

			}
		}

		fclose(fp_inicialization_data);
	}
	printf("save_velocity_for_init Ok\n");
} // save_velocity_for_init


doublereal signlog10(doublereal r21) {
	if (r21 > 2.0) {
		return log10(r21);
	}
	else if (r21 < -2.0) {
		return -log10(fabs(r21));
	}
	else {
		return 0.0;
	}
}

// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2_apparat_hot( integer maxelm, integer ncell,
	FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate,
	bool bextendedprint, integer ikey, BLOCK* &b)
{
	const bool lite_export = true;
	// 16 знаков после запятой сохранять никому ненужно,
	// вполне достаточно шести знаков.

#if doubleintprecision == 1
	printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
#else
	printf("ionly_solid_visible =%d\n", ionly_solid_visible);
#endif
	

	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal* temp_shadow = nullptr;
	if (ionly_solid_visible == 1) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == 1) {
		total_deformation_shadow = new doublereal*[4];
		for (integer j_1 = 0; j_1 < 4; j_1++) {
			total_deformation_shadow[j_1] = new doublereal[t.maxelm + t.maxbound];
			for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
				total_deformation_shadow[j_1][i_1] = t.total_deformation[j_1][i_1];
			}
		}
	}


	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	err = 0;
	switch (ikey) {
		case 0:  fp=fopen64("ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
		case 1:  fp=fopen64("ALICEFLOW0_27_temp_apparat_hot.PLT", "wb");  break; // То что нужно для отчёта Алексею.
		default: fp=fopen64("ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
    }
	if (fp == NULL) err = 1;
#else
	switch (ikey) {
	case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
	case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_apparat_hot.PLT", "wb");  break; // То что нужно для отчёта Алексею.
	default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
	}
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in exporttecplotxy360T_3D_part2_apparat_hot in my_export_tecplot.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1=fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1 == nullptr) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		if (bOk)
		{


			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // считаем что температура всегда активна
				}
			}

			if (ivarexport == 1) {

				integer ncell_shadow = 0;
				for (i = 0; i < ncell; i++) {

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					//in_model_temp(p1, ib1, b, lb);
					//in_model_temp(p2, ib2, b, lb);
					//in_model_temp(p3, ib3, b, lb);
					//in_model_temp(p4, ib4, b, lb);
					//in_model_temp(p5, ib5, b, lb);
					//in_model_temp(p6, ib6, b, lb);
					//in_model_temp(p7, ib7, b, lb);
					//in_model_temp(p8, ib8, b, lb);

					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];

					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
						ncell_shadow++;
					}
				}
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// запись имён переменных
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
				fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // печатается только поле температур

#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif

				

				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			else if (ivarexport == 3) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ничего не делаем
					ncell_shadow = 0;
					for (i = 0; i < ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						/*
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);
						*/

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
							ncell_shadow++;
						}
					}
					//ncell_shadow = ncell;
				}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						ncell_shadow = 0;
						for (i = 0; i < t.database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = t.database.nvtxcell[0][i] - 1;
							inode2 = t.database.nvtxcell[1][i] - 1;
							inode3 = t.database.nvtxcell[2][i] - 1;
							inode4 = t.database.nvtxcell[3][i] - 1;
							inode5 = t.database.nvtxcell[4][i] - 1;
							inode6 = t.database.nvtxcell[5][i] - 1;
							inode7 = t.database.nvtxcell[6][i] - 1;
							inode8 = t.database.nvtxcell[7][i] - 1;
							
							integer inode2W= t.sosedi[WSIDE][inode1].iNODE1;
							integer inode3W= t.sosedi[WSIDE][inode4].iNODE1;
							integer inode6W= t.sosedi[WSIDE][inode5].iNODE1;
							integer inode7W= t.sosedi[WSIDE][inode8].iNODE1;
							
							
							integer inode5B= t.sosedi[BSIDE][inode1].iNODE1;
							integer inode6B= t.sosedi[BSIDE][inode2].iNODE1;
							integer inode7B= t.sosedi[BSIDE][inode3].iNODE1;
							integer inode8B= t.sosedi[BSIDE][inode4].iNODE1;
							

							
							integer inode3S = t.sosedi[SSIDE][inode2].iNODE1; 
							integer inode4S = t.sosedi[SSIDE][inode1].iNODE1; 
							integer inode7S = t.sosedi[SSIDE][inode6].iNODE1; 
							integer inode8S = t.sosedi[SSIDE][inode5].iNODE1; 
							
							//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
							//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
							//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
							//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
							//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
							//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
							//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
							//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
							//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
							//pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
							//pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
							//pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
							//in_model_temp(p1, ib1, b, lb);
							//in_model_temp(p2, ib2, b, lb);
							//in_model_temp(p3, ib3, b, lb);
							//in_model_temp(p4, ib4, b, lb);
							//in_model_temp(p5, ib5, b, lb);
							//in_model_temp(p6, ib6, b, lb);
							//in_model_temp(p7, ib7, b, lb);
							//in_model_temp(p8, ib8, b, lb);

							ib1 = t.whot_is_block[inode1];
							ib2 = t.whot_is_block[inode2];
							ib3 = t.whot_is_block[inode3];
							ib4 = t.whot_is_block[inode4];
							ib5 = t.whot_is_block[inode5];
							ib6 = t.whot_is_block[inode6];
							ib7 = t.whot_is_block[inode7];
							ib8 = t.whot_is_block[inode8];


							// визуализация только твёрдого тела.
							// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
							if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
								&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1))) 
							{
								if ((b[ib1].bvisible) && (b[ib2].bvisible)&& (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
									ncell_shadow++;
								}
								
							}
							else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
							{
								
								if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
									ncell_shadow++;
									temp_shadow[inode5] = t.potent[inode1];
									temp_shadow[inode6] = t.potent[inode2];
									temp_shadow[inode7] = t.potent[inode3];
									temp_shadow[inode8] = t.potent[inode4];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[j_4][inode5] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode2];
										total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode3];
										total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode4];
									}
								}
							}
							else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
							{
								
								if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

									ncell_shadow++;
									temp_shadow[inode2] = t.potent[inode1];
									temp_shadow[inode3] = t.potent[inode4];
									temp_shadow[inode6] = t.potent[inode5];
									temp_shadow[inode7] = t.potent[inode8];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[j_4][inode2] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode4];
										total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode5];
										total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode8];
									}
								}
							}
							else if (((inode3S >= 0) && (inode3S < t.maxelm)&& (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
                            {
								
								if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

									ncell_shadow++;
									temp_shadow[inode4] = t.potent[inode1];
									temp_shadow[inode3] = t.potent[inode2];
									temp_shadow[inode8] = t.potent[inode5];
									temp_shadow[inode7] = t.potent[inode6];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[j_4][inode4] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode2];
										total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode5];
										total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode6];
									}
								}
							}
						}
						if (ncell_shadow == 0) {
							// У нас совсем нету твёрдого тела, поэтому мы показываем только  жидкость.
							ncell_shadow = ncell;
							ionly_solid_visible = 0;
						}
					}
				}
				// Полный набор искомых величин и теплопередача и гидродинамика:
				//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
				fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				

#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif

				
				

				if (bvery_big_memory) {
					if (lite_export) {
						// запись x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					else {
						// запись x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части

		// Запись поля температур производится всегда.


		// запись температуры
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == 1) {
					fprintf(fp, "%+.3f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.3f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.3f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.3f ", t.potent[i]);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}
		else {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == 1) {
					fprintf(fp, "%+.16f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.16f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.16f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.16f ", t.potent[i]);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		


		fprintf(fp, "\n");

		// Lam
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.4f ", t.prop[LAM][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < t.maxbound; i++) {
					fprintf(fp, "%+.4f ", t.prop_b[LAM][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}
		else {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", t.prop[LAM][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < t.maxbound; i++) {
					fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		fprintf(fp, "\n");

		// Запись гидродинамических величин если необходимо:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VX][i + maxelm] * f[idfluid].potent[VX][i + maxelm];
					doublereal svy = f[idfluid].potent[VY][i + maxelm] * f[idfluid].potent[VY][i + maxelm];
					doublereal svz = f[idfluid].potent[VZ][i + maxelm] * f[idfluid].potent[VZ][i + maxelm];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}


			fprintf(fp, "\n");

			// Pressure
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Pressure
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// PAM
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// PAM
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				}
				else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Rho
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mu
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // Турбулентная динамическая вязкость
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]/ f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					
					//fprintf(fp, "%+.16f ", 1.0*f[t.ptr[1][i]].icolor_different_fluid_domain[t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// MUT
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					//fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]/ f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime == ZEROEQMOD) || 
						(f[0].iflowregime == SMAGORINSKY)||
						(f[0].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == RANS_MENTER_SST) ||
						(f[0].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // Завихрённость - модуль ротора скорости.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Curl
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i + maxelm]); // GRADZVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

		}

			doublereal *Tx = nullptr;
			doublereal *Ty = nullptr;
			doublereal *Tz = nullptr;
			Tx = new doublereal[t.maxelm + t.maxbound];
			Ty = new doublereal[t.maxelm + t.maxbound];
			Tz = new doublereal[t.maxelm + t.maxbound];

			// инициализация нулём.
			for (i = 0; i<t.maxelm + t.maxbound; i++) {
				Tx[i] = 0.0;
				Ty[i] = 0.0;
				Tz[i] = 0.0;
			}

			// нахождение градиентов.
			for (i = 0; i<t.maxelm; i++) {
				// Только внутренние узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, false,
					t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
			}

			for (i = 0; i<t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, true,
					t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
			}


			// Сохранение в файл.

			if (lite_export) {
				doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

				// Heat Flux X
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1+1;
						buf1 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1+2;
						buf2 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 3;
						buf3 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 4;
						buf4 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 5;
						buf5 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 6;
						buf6 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 7;
						buf7 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 8;
						buf8 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						i = i1 + 9;
						buf9 = signlog10 (-t.prop[LAM][i] * Tx[i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10 (-t.prop[LAM][i] * Tx[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+1;
							buf1 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+2;
							buf2 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+3;
							buf3 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+4;
							buf4 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+5;
							buf5 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+6;
							buf6 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+7;
							buf7 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+8;
							buf8 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1+9;
							buf9 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Y
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop[LAM][i] * Ty[i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Z
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop[LAM][i] * Tz[i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// Mag Heat Flux
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1+9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							 i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1+7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1+9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0,buf1,buf2,buf3,buf4,buf5,buf6,buf7,buf8,buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// log10 Mag Heat Flux
				const doublereal eps_min = 2.0;
				const doublereal d_stub = 0.0;
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}
							if (buf1 > eps_min) {
								buf1 = log10(buf1);
							}
							else {
								buf1 = d_stub;
							}
							if (buf2 > eps_min) {
								buf2 = log10(buf2);
							}
							else {
								buf2 = d_stub;
							}
							if (buf3 > eps_min) {
								buf3 = log10(buf3);
							}
							else {
								buf3 = d_stub;
							}
							if (buf4 > eps_min) {
								buf4 = log10(buf4);
							}
							else {
								buf4 = d_stub;
							}
							if (buf5 > eps_min) {
								buf5 = log10(buf5);
							}
							else {
								buf5 = d_stub;
							}
							if (buf6 > eps_min) {
								buf6 = log10(buf6);
							}
							else {
								buf6 = d_stub;
							}
							if (buf7 > eps_min) {
								buf7 = log10(buf7);
							}
							else {
								buf7 = d_stub;
							}
							if (buf8 > eps_min) {
								buf8 = log10(buf8);
							}
							else {
								buf8 = d_stub;
							}
							if (buf9 > eps_min) {
								buf9 = log10(buf9);
							}
							else {
								buf9 = d_stub;
							}

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}

							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

			}
			else {
				doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

				// Heat Flux X
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 1;
						buf1 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 2;
						buf2 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 3;
						buf3 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 4;
						buf4 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 5;
						buf5 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 6;
						buf6 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 7;
						buf7 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 8;
						buf8 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 9;
						buf9 = -t.prop[LAM][i] * Tx[i];

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop[LAM][i] * Tx[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 1;
							buf1 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 2;
							buf2 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 3;
							buf3 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 4;
							buf4 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 5;
							buf5 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 6;
							buf6 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 7;
							buf7 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 8;
							buf8 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 9;
							buf9 = -t.prop_b[LAM][i] * Tx[i + maxelm];


							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							fprintf(fp, "%+.16f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Y
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 1;
						buf1 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 2;
						buf2 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 3;
						buf3 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 4;
						buf4 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 5;
						buf5 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 6;
						buf6 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 7;
						buf7 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 8;
						buf8 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 9;
						buf9 = -t.prop[LAM][i] * Ty[i];

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop[LAM][i] * Ty[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 1;
							buf1 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 2;
							buf2 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 3;
							buf3 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 4;
							buf4 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 5;
							buf5 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 6;
							buf6 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 7;
							buf7 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 8;
							buf8 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 9;
							buf9 = -t.prop_b[LAM][i] * Ty[i + maxelm];


							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							fprintf(fp, "%+.16f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Z
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 1;
						buf1 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 2;
						buf2 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 3;
						buf3 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 4;
						buf4 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 5;
						buf5 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 6;
						buf6 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 7;
						buf7 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 8;
						buf8 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 9;
						buf9 = -t.prop[LAM][i] * Tz[i];

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop[LAM][i] * Tz[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 1;
							buf1 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 2;
							buf2 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 3;
							buf3 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 4;
							buf4 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 5;
							buf5 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 6;
							buf6 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 7;
							buf7 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 8;
							buf8 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 9;
							buf9 = -t.prop_b[LAM][i] * Tz[i + maxelm];


							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							fprintf(fp, "%+.16f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// Mag Heat Flux
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							 i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							 i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.16f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");
				
				
				// log10 Mag Heat Flux
				const doublereal eps_min = 2.0;
				const doublereal d_stub = 0.0;
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}
							if (buf1 > eps_min) {
								buf1 = log10(buf1);
							}
							else {
								buf1 = d_stub;
							}
							if (buf2 > eps_min) {
								buf2 = log10(buf2);
							}
							else {
								buf2 = d_stub;
							}
							if (buf3 > eps_min) {
								buf3 = log10(buf3);
							}
							else {
								buf3 = d_stub;
							}
							if (buf4 > eps_min) {
								buf4 = log10(buf4);
							}
							else {
								buf4 = d_stub;
							}
							if (buf5 > eps_min) {
								buf5 = log10(buf5);
							}
							else {
								buf5 = d_stub;
							}
							if (buf6 > eps_min) {
								buf6 = log10(buf6);
							}
							else {
								buf6 = d_stub;
							}
							if (buf7 > eps_min) {
								buf7 = log10(buf7);
							}
							else {
								buf7 = d_stub;
							}
							if (buf8 > eps_min) {
								buf8 = log10(buf8);
							}
							else {
								buf8 = d_stub;
							}
							if (buf9 > eps_min) {
								buf9 = log10(buf9);
							}
							else {
								buf9 = d_stub;
							}

							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							
							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}

							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.16f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			
				fprintf(fp, "\n");


			}
			

			// Освобождение оперативной памяти.
			if (Tx != nullptr) {
				delete[] Tx;
			}
			if (Ty != nullptr) {
				delete[] Ty;
			}
			if (Tz != nullptr) {
				delete[] Tz;
			}

			fprintf(fp, "\n");

			for (integer j_6 = 0; j_6 < 4; j_6++) {

				// запись полной деформации
				if (lite_export) {
					for (i = 0; i < maxelm; i++) {


						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							//printf("%e \n", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							//printf("%e \n", t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) {
							fprintf(fp, "\n");
							//getchar();
						}
					}

					if (bextendedprint) {
						for (i = maxelm; i < maxelm + t.maxbound; i++) {
							if (ionly_solid_visible == 1) {
								fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							}
							else {
								fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							}
							if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					for (i = 0; i < maxelm; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (bextendedprint) {
						for (i = maxelm; i < maxelm + t.maxbound; i++) {
							if (ionly_solid_visible == 1) {
								fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							}
							else {
								fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							}
							if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}


				fprintf(fp, "\n");

			}

		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					/*
					TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
					*/
					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					/*
					in_model_temp(p1, ib1, b, lb);
					in_model_temp(p2, ib2, b, lb);
					in_model_temp(p3, ib3, b, lb);
					in_model_temp(p4, ib4, b, lb);
					in_model_temp(p5, ib5, b, lb);
					in_model_temp(p6, ib6, b, lb);
					in_model_temp(p7, ib7, b, lb);
					in_model_temp(p8, ib8, b, lb);
					*/

					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];


					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						
					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
						integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
						integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
						integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


						integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
						integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
						integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
						integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



						integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
						integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
						integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
						integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;

						
						//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
						//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						//pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
						//pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
						//pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						//in_model_temp(p1, ib1, b, lb);
						//in_model_temp(p2, ib2, b, lb);
						//in_model_temp(p3, ib3, b, lb);
						//in_model_temp(p4, ib4, b, lb);
						//in_model_temp(p5, ib5, b, lb);
						//in_model_temp(p6, ib6, b, lb);
						//in_model_temp(p7, ib7, b, lb);
						//in_model_temp(p8, ib8, b, lb);

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];

						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{ 
							
							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

									}
						}
						else if ( ((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

								}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

									}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
					    {
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

								}
						}
					}
					else {
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						//in_model_temp(p1, ib1, b, lb);
						//in_model_temp(p2, ib2, b, lb);
						//in_model_temp(p3, ib3, b, lb);
						//in_model_temp(p4, ib4, b, lb);
						//in_model_temp(p5, ib5, b, lb);
						//in_model_temp(p6, ib6, b, lb);
						//in_model_temp(p7, ib7, b, lb);
						//in_model_temp(p8, ib8, b, lb);

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							
						}

					}
				}
			}
		}
		else {
#ifdef MINGW_COMPILLER
		err = 0;
		fp1 = fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
		if (fp1 == nullptr) err = 1;
#else
		err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	if (temp_shadow != nullptr) {
		delete[] temp_shadow;
		temp_shadow = nullptr;
	}
	//total_deformation
	
	if (total_deformation_shadow != nullptr) {
		for (integer j_6 = 0; j_6 < 4; j_6++) {
			if (total_deformation_shadow[j_6] != nullptr) {
				delete[] total_deformation_shadow[j_6];
				total_deformation_shadow[j_6] = nullptr;
			}
		}
		delete[] total_deformation_shadow;
		total_deformation_shadow = nullptr;
	}
	

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);



}

// 4.01.2018
// Запрограммировано за один вечер.
// Данная функция необходима для подготовки картинок для отчетов на основе огромных моделей.
// Т.к. её аналог на делпфи совершенно не справляется с большими файлами по причине ограниченных функциональных возможностей библиоек стандартных компонентов на очень больших размерах.
void tecplot360patcher_for_print_in_report() {
	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	if (1 && steady_or_unsteady_global_determinant == 5) {
		// Второй температурный солвер.
		fp = fopen64("ALICEFLOW0_08_temp.PLT", "r");
    }
	else {
		fp = fopen64("ALICEFLOW0_07_temp.PLT", "r");
	}
	if (fp==NULL) err = 1;
#else
	if (1 && steady_or_unsteady_global_determinant == 5) {
		// Второй температурный солвер.
		err = fopen_s(&fp, "ALICEFLOW0_08_temp.PLT", "r");
}
	else {
		err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "r");
	}
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in function tecplot360patcher_for_print_in_report in my_export_tecplotr3.c\n");
		//getchar();
		system("pause");

	}
	else {
		FILE *fp1=NULL;
		errno_t err1=0;
#ifdef MINGW_COMPILLER
		fp1 = fopen64("ALICEFlow0_07_Visualisation_Magement", "w");
		if (fp1==NULL) err1 = 1;
#else
		err1 = fopen_s(&fp1, "ALICEFlow0_07_Visualisation_Magement", "w");
#endif
		
		if ((err1) != 0) {
			printf("Create File ALICEFlow0_07_Visualisation_Magement Error\n");
			//getchar();
			system("pause");

		}
		else {

			integer iN = -1, iE = -1, iVar = -1;
			char str[600] = "\0";

			//fscanf_s(fp, "%s", str);
			fgets(str, sizeof(str), fp);
			printf("%s", str);
			fprintf(fp1, "%s", str);
			//str = "\0";
			//fscanf_s(fp, "%s", str);
			fgets(str, sizeof(str), fp);
			printf("%s", str);
			fprintf(fp1, "%s", str);
			while (strlen(str) == 0) {
				//str = "\0";
				//fscanf_s(fp, "%s", str);
				fgets(str, sizeof(str), fp);
				printf("%s", str);
				fprintf(fp1, "%s", str);
			}
			if (strlen(str) > 0) {
				//printf("%s \n",str);
				//printf("strlen=%lld\n", strlen(str));
				while (strlen(str) == 1) {
					fgets(str, sizeof(str), fp);
					printf("%s", str);
					fprintf(fp1, "%s", str);
				}
				iVar++;
				//printf("iVar=%lld\n", iVar);
				//system("pause");
				size_t length_str = strlen(str);
				for (int i = 0; i < length_str; i++) if (str[i] == ',') {
					iVar++;
					//printf("i=%d iVar=%lld\n",i, iVar);
					//system("pause");
				}
				iVar++;
				//printf("iVar=%lld\n", iVar);
				//system("pause");
			}
			//str = "\0";
			//fscanf_s(fp, "%s", str);
			fgets(str, sizeof(str), fp);
			printf("%s", str);
			while (strlen(str) == 0) {
				fgets(str, sizeof(str), fp);
			}
			if (strlen(str) > 0) {
				integer ic1 = 0;
				bool bf11 = false;
				bool bf1 = true;
				size_t length_str = strlen(str);
				for (int i = 0; i < length_str; i++) if (str[i] == '=') {

					if (ic1 == 1 || ic1 == 2) {
						bf11 = true;
					}
					else {
						bf11 = false;
					}
					if (bf11&&bf1) {
						bf1 = false;
						char str1[600] = "\0";
						integer k = 0;
						//size_t length_str = strlen(str);
						for (int j = i + 1; j < length_str; j++) {
							if ((str[j] >= '0') && (str[j] <= '9')) {
								str1[k++] = str[j];
							}
							else {
								str1[k] = '\0';
								iN = atol(str1);
								//printf("%s %lld\n",str1,iN);
								break;
							}
						}
					}
					else if (bf11 && !bf1) {
						char str1[600] = "\0";
						integer k = 0;
						//size_t length_str = strlen(str);
						for (int j = i + 1; j < length_str; j++) {
							if ((str[j] >= '0') && (str[j] <= '9')) {
								str1[k++] = str[j];
							}
							else {
								str1[k] = '\0';
								iE = atol(str1);
								//printf("%s %lld\n", str1, iE);
								break;
							}
						}
					}
					ic1++;
				}
			}
			printf("iN=%lld iE=%lld number Variable =%lld\n", iN, iE, iVar);
			//system("pause");

			doublereal** func = nullptr;
			func = new doublereal*[iVar];
			for (integer i = 0; i < iVar; i++) func[i] = new doublereal[iN];

			// Считывание всех функций.
			for (integer i = 0; i < iVar; i++) {
				for (integer j = 0; j < iN; j++) {
					float fin = 0.0;
#ifdef MINGW_COMPILLER
					fscanf(fp, "%f", &fin);
#else
					fscanf_s(fp, "%f", &fin);
#endif
					
					func[i][j] = fin;
				}
				printf("compleate %lld is %lld\n", i + 1, iVar);
			}

			//printf("1 - X, 2 - Y, 3 - Z\n");
			//char ch2 = getchar();
			char ch2 = '2'; // инициализируем осью Y
			switch (pfpir.idir) {
			case 0: ch2 = '1'; break;
			case 1: ch2 = '2'; break;
			case 2: ch2 = '3'; break;
			default :
				printf("Unknown directional in variable pfpir.idir=%lld\n", pfpir.idir);
				printf("1 - X, 2 - Y, 3 - Z - anather unknown.\n");
				printf("ERROR function tecplot360pather_for_print_in_report();");
				system("PAUSE");
				exit(1);
				break;
			}
			//printf("print minimum = \n");
			doublereal fmin1 = 0.0;
			//scanf("%f", &fmin1);
			fmin1 = pfpir.fminimum;
			//printf("print maximum = \n");
			doublereal fmax1 = 0.0;
			//scanf("%f", &fmax1);
			fmax1 = pfpir.fmaximum;

			switch (ch2) {
			case '1': printf("X directional interval minimum =%e maximum %e\n", fmin1, fmax1);
				break;
			case '2': printf("Y directional interval minimum =%e maximum %e\n", fmin1, fmax1);
				break;
			case '3': printf("Z directional interval minimum =%e maximum %e\n", fmin1, fmax1);
				break;
			}
			//system("pause");

			integer** nvtx_1 = nullptr;
			nvtx_1 = new integer*[iE];
			for (integer i = 0; i < iE; i++) nvtx_1[i] = new integer[8];

			integer iE2 = 0;
			for (integer i = 0; i < iE; i++) {
				integer ie1 = -1;
				integer ie2 = -1;
				integer ie3 = -1;
				integer ie4 = -1;
				integer ie5 = -1;
				integer ie6 = -1;
				integer ie7 = -1;
				integer ie8 = -1;
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
				fscanf(fp, "%lld", &ie1);
				ie1--;
				fscanf(fp, "%lld", &ie2);
				ie2--;
				fscanf(fp, "%lld", &ie3);
				ie3--;
				fscanf(fp, "%lld", &ie4);
				ie4--;
				fscanf(fp, "%lld", &ie5);
				ie5--;
				fscanf(fp, "%lld", &ie6);
				ie6--;
				fscanf(fp, "%lld", &ie7);
				ie7--;
				fscanf(fp, "%lld", &ie8);
				ie8--;
#else
				fscanf(fp, "%d", &ie1);
				ie1--;
				fscanf(fp, "%d", &ie2);
				ie2--;
				fscanf(fp, "%d", &ie3);
				ie3--;
				fscanf(fp, "%d", &ie4);
				ie4--;
				fscanf(fp, "%d", &ie5);
				ie5--;
				fscanf(fp, "%d", &ie6);
				ie6--;
				fscanf(fp, "%d", &ie7);
				ie7--;
				fscanf(fp, "%d", &ie8);
				ie8--;
#endif
#else
#if doubleintprecision == 1
				fscanf_s(fp, "%lld", &ie1);
				ie1--;
				fscanf_s(fp, "%lld", &ie2);
				ie2--;
				fscanf_s(fp, "%lld", &ie3);
				ie3--;
				fscanf_s(fp, "%lld", &ie4);
				ie4--;
				fscanf_s(fp, "%lld", &ie5);
				ie5--;
				fscanf_s(fp, "%lld", &ie6);
				ie6--;
				fscanf_s(fp, "%lld", &ie7);
				ie7--;
				fscanf_s(fp, "%lld", &ie8);
				ie8--;
#else
				fscanf_s(fp, "%d", &ie1);
				ie1--;
				fscanf_s(fp, "%d", &ie2);
				ie2--;
				fscanf_s(fp, "%d", &ie3);
				ie3--;
				fscanf_s(fp, "%d", &ie4);
				ie4--;
				fscanf_s(fp, "%d", &ie5);
				ie5--;
				fscanf_s(fp, "%d", &ie6);
				ie6--;
				fscanf_s(fp, "%d", &ie7);
				ie7--;
				fscanf_s(fp, "%d", &ie8);
				ie8--;
#endif
#endif

				nvtx_1[i][0] = ie1 + 1;
				nvtx_1[i][1] = ie2 + 1;
				nvtx_1[i][2] = ie3 + 1;
				nvtx_1[i][3] = ie4 + 1;
				nvtx_1[i][4] = ie5 + 1;
				nvtx_1[i][5] = ie6 + 1;
				nvtx_1[i][6] = ie7 + 1;
				nvtx_1[i][7] = ie8 + 1;
				doublereal xo = 0.0;
				doublereal yo = 0.0;
				doublereal zo = 0.0;
				xo = 0.125*(func[0][ie1] + func[0][ie2] + func[0][ie3] + func[0][ie4] + func[0][ie5] + func[0][ie6] + func[0][ie7] + func[0][ie8]);
				yo = 0.125*(func[1][ie1] + func[1][ie2] + func[1][ie3] + func[1][ie4] + func[1][ie5] + func[1][ie6] + func[1][ie7] + func[1][ie8]);
				zo = 0.125*(func[2][ie1] + func[2][ie2] + func[2][ie3] + func[2][ie4] + func[2][ie5] + func[2][ie6] + func[2][ie7] + func[2][ie8]);
				switch (ch2) {
				case '1':if ((xo >= fmin1) && (xo <= fmax1)) iE2++;
					break;
				case '2':if ((yo >= fmin1) && (yo <= fmax1)) iE2++;
					break;
				case '3':if ((zo >= fmin1) && (zo <= fmax1)) iE2++;
					break;
				}
			}

#if doubleintprecision == 1
			fprintf(fp1, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", iN, iE2);
#else
			fprintf(fp1, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", iN, iE2);
#endif


			for (integer i = 0; i < iVar; i++) {
				for (integer j = 0; j < iN; j++) {

					fprintf(fp1, "%e ", func[i][j]);
					if (j % 20 == 0) fprintf(fp1, "\n");
				}
				fprintf(fp1, "\n");
				printf("compleate %lld is %lld\n", i + 1, iVar);
			}

			for (integer i = 0; i < iE; i++) {
				integer ie1 = nvtx_1[i][0];
				integer ie2 = nvtx_1[i][1];
				integer ie3 = nvtx_1[i][2];
				integer ie4 = nvtx_1[i][3];
				integer ie5 = nvtx_1[i][4];
				integer ie6 = nvtx_1[i][5];
				integer ie7 = nvtx_1[i][6];
				integer ie8 = nvtx_1[i][7];
				ie1--;
				ie2--;
				ie3--;
				ie4--;
				ie5--;
				ie6--;
				ie7--;
				ie8--;
				doublereal xo = 0.0;
				doublereal yo = 0.0;
				doublereal zo = 0.0;
				xo = 0.125*(func[0][ie1] + func[0][ie2] + func[0][ie3] + func[0][ie4] + func[0][ie5] + func[0][ie6] + func[0][ie7] + func[0][ie8]);
				yo = 0.125*(func[1][ie1] + func[1][ie2] + func[1][ie3] + func[1][ie4] + func[1][ie5] + func[1][ie6] + func[1][ie7] + func[1][ie8]);
				zo = 0.125*(func[2][ie1] + func[2][ie2] + func[2][ie3] + func[2][ie4] + func[2][ie5] + func[2][ie6] + func[2][ie7] + func[2][ie8]);

#if doubleintprecision == 1
				switch (ch2) {
				case '1':if ((xo >= fmin1) && (xo <= fmax1)) {
					fprintf(fp1, "%lld %lld %lld %lld %lld %lld %lld %lld\n", ie1 + 1, ie2 + 1, ie3 + 1, ie4 + 1, ie5 + 1, ie6 + 1, ie7 + 1, ie8 + 1);
				}
						 break;
				case '2':if ((yo >= fmin1) && (yo <= fmax1)) {
					fprintf(fp1, "%lld %lld %lld %lld %lld %lld %lld %lld\n", ie1 + 1, ie2 + 1, ie3 + 1, ie4 + 1, ie5 + 1, ie6 + 1, ie7 + 1, ie8 + 1);
				}
						 break;
				case '3':if ((zo >= fmin1) && (zo <= fmax1)) {
					fprintf(fp1, "%lld %lld %lld %lld %lld %lld %lld %lld\n", ie1 + 1, ie2 + 1, ie3 + 1, ie4 + 1, ie5 + 1, ie6 + 1, ie7 + 1, ie8 + 1);
				}
						 break;
				}
#else
				switch (ch2) {
				case '1':if ((xo >= fmin1) && (xo <= fmax1)) {
					fprintf(fp1, "%d %d %d %d %d %d %d %d\n", ie1 + 1, ie2 + 1, ie3 + 1, ie4 + 1, ie5 + 1, ie6 + 1, ie7 + 1, ie8 + 1);
				}
						 break;
				case '2':if ((yo >= fmin1) && (yo <= fmax1)) {
					fprintf(fp1, "%d %d %d %d %d %d %d %d\n", ie1 + 1, ie2 + 1, ie3 + 1, ie4 + 1, ie5 + 1, ie6 + 1, ie7 + 1, ie8 + 1);
				}
						 break;
				case '3':if ((zo >= fmin1) && (zo <= fmax1)) {
					fprintf(fp1, "%d %d %d %d %d %d %d %d\n", ie1 + 1, ie2 + 1, ie3 + 1, ie4 + 1, ie5 + 1, ie6 + 1, ie7 + 1, ie8 + 1);
				}
						 break;
			}
#endif
			}
			fprintf(fp1, "\n");

			// Освобождение оперативной памяти.
			for (integer i = 0; i < iE; i++) {
				delete[] nvtx_1[i];
				nvtx_1[i] = nullptr;
			}
			delete[] nvtx_1;
			nvtx_1 = nullptr;


			for (integer i = 0; i < iVar; i++) {
				delete[] func[i]; func[i] = nullptr;
		    }

			delete[] func;
			func = nullptr;

			fclose(fp);
			fclose(fp1);

			if (1 && steady_or_unsteady_global_determinant == 5) {
				printf("file ALICEFLOW0_08_temp.PLT preobrazovan. Ok.\n");
			}
			else {
				printf("file ALICEFLOW0_07_temp.PLT preobrazovan. Ok.\n");
			}
			//system("pause");

		}
	}

}// tecplot360patcher_for_print_in_report


// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t,
	integer flow_interior_count, integer ianimate, bool bextendedprint,
	integer ikey, BLOCK* &b, integer &lb)
{

	bool bFLAGitoa = false;

	const bool lite_export = true;
	// 16 знаков после запятой сохранять никому ненужно,
	// вполне достаточно шести знаков.

#if doubleintprecision == 1
	printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
#else
	printf("ionly_solid_visible =%d\n", ionly_solid_visible);
#endif


	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp=NULL; // Выходной файл.
	FILE *fp1=NULL; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal* temp_shadow = nullptr;
	if (ionly_solid_visible == 1) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == 1) {
		total_deformation_shadow = new doublereal*[4];
		for (integer j_1 = 0; j_1 < 4; j_1++) {
			total_deformation_shadow[j_1] = new doublereal[t.maxelm + t.maxbound];
			for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
				total_deformation_shadow[j_1][i_1] = t.total_deformation[j_1][i_1];
			}
		}
	}


	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	err = 0;
	switch (ikey) {
	case 0:  fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1:  fp = fopen64("ALICEFLOW0_27_temp.PLT", "wb");  break; // То что нужно для отчёта Алексею.
	default: fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
}
	if (fp==NULL) err = 1;
#else
	switch (ikey) {
	case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp.PLT", "wb");  break; // То что нужно для отчёта Алексею.
	default: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	}
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2 in my_export_tecplot3.c\n");
		system("pause");

	}
	else {

		int c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1=fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1==NULL) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		if (bOk)
		{


			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // считаем что температура всегда активна
				}
			}

			if (ivarexport == 1) {

				integer ncell_shadow = 0;
				for (i = 0; i < ncell; i++) {

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					// Использование быстродействующей хеш таблицы whot_is_block значительно
					// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
					ib1 = t.whot_is_block[inode1];
					//in_model_temp(p1, ib1, b, lb);	
					//if (ib1 != t.whot_is_block[inode1]) {
						//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
						//getchar();
					//}
					ib2 = t.whot_is_block[inode2];
					//in_model_temp(p2, ib2, b, lb);
					//if (ib2 != t.whot_is_block[inode2]) {
						//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
						//getchar();
					//}
					ib3 = t.whot_is_block[inode3];
					//in_model_temp(p3, ib3, b, lb);
					//if (ib3 != t.whot_is_block[inode3]) {
						//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
						//getchar();
					//}
					ib4 = t.whot_is_block[inode4];
					//in_model_temp(p4, ib4, b, lb);
					//if (ib4 != t.whot_is_block[inode4]) {
						//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
						//getchar();
					//}
					ib5 = t.whot_is_block[inode5];
					//in_model_temp(p5, ib5, b, lb);
					//if (ib5 != t.whot_is_block[inode5]) {
						//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
						//getchar();
					//}
					ib6 = t.whot_is_block[inode6];
					//in_model_temp(p6, ib6, b, lb);
					//if (ib6 != t.whot_is_block[inode6]) {
						//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
						//getchar();
					//}
					ib7 = t.whot_is_block[inode7];
					//in_model_temp(p7, ib7, b, lb);
					//if (ib7 != t.whot_is_block[inode7]) {
						//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
						//getchar();
					//}
					ib8 = t.whot_is_block[inode8];
					//in_model_temp(p8, ib8, b, lb);
					//if (ib8 != t.whot_is_block[inode8]) {
						//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
						//getchar();
					//}
					

					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
						ncell_shadow++;
					}
				}
				
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// запись имён переменных
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
				if (1 && 2 == steady_or_unsteady_global_determinant) {
					// Вызов только сеточного генератора
					fprintf(fp, "VARIABLES = x, y, z, quolity, log(quolity) \n");  // печатается только расчётная сетка
				}
				else {
					// печатается только поле температур
					fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // печатается только поле температур
				}

#if doubleintprecision == 1
																																																							   // запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
																																																							   // запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif



				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (1 && 2 == steady_or_unsteady_global_determinant) {
						// quolity
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
							volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", dmax / dmin);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						// log10(quolity)
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
							volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", log10(dmax / dmin));
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			else if (ivarexport == 3) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ничего не делаем
					ncell_shadow = 0;
					for (i = 0; i < ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						/*
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);
						*/

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) &&
							(b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
							ncell_shadow++;
						}
					}
					//ncell_shadow = ncell;
					
				}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						ncell_shadow = 0;
						for (i = 0; i < t.database.ncell; i++) {

							

								integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
								inode1 = t.database.nvtxcell[0][i] - 1;
								inode2 = t.database.nvtxcell[1][i] - 1;
								inode3 = t.database.nvtxcell[2][i] - 1;
								inode4 = t.database.nvtxcell[3][i] - 1;
								inode5 = t.database.nvtxcell[4][i] - 1;
								inode6 = t.database.nvtxcell[5][i] - 1;
								inode7 = t.database.nvtxcell[6][i] - 1;
								inode8 = t.database.nvtxcell[7][i] - 1;

								if (b[t.whot_is_block[inode1]].bvisible&&
									b[t.whot_is_block[inode2]].bvisible&&
									b[t.whot_is_block[inode3]].bvisible&&
									b[t.whot_is_block[inode4]].bvisible&&
									b[t.whot_is_block[inode5]].bvisible&&
									b[t.whot_is_block[inode6]].bvisible&&
									b[t.whot_is_block[inode7]].bvisible&&
									b[t.whot_is_block[inode8]].bvisible) {

								integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
								integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
								integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
								integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


								integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
								integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
								integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
								integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



								integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
								integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
								integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
								integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;

								TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
								center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
								center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
								center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
								center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
								center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
								center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
								center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
								center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
								pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
								pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
								pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

								integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
								in_model_temp(p1, ib1, b, lb);
								in_model_temp(p2, ib2, b, lb);
								in_model_temp(p3, ib3, b, lb);
								in_model_temp(p4, ib4, b, lb);
								in_model_temp(p5, ib5, b, lb);
								in_model_temp(p6, ib6, b, lb);
								in_model_temp(p7, ib7, b, lb);
								in_model_temp(p8, ib8, b, lb);

								// визуализация только твёрдого тела.
								// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
								if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
									&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
								{
									if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) &&
										(b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
										ncell_shadow++;
									}

								}
								else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) &&
									(inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) &&
									(t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) &&
									(t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) &&
									(t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) &&
										(t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
								{

									//if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible))
									if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) &&
										(b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible))
									{
										ncell_shadow++;
										temp_shadow[inode5] = t.potent[inode1];
										temp_shadow[inode6] = t.potent[inode2];
										temp_shadow[inode7] = t.potent[inode3];
										temp_shadow[inode8] = t.potent[inode4];
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[j_4][inode5] = t.total_deformation[j_4][inode1];
											total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode2];
											total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode3];
											total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode4];
										}
									}
								}
								else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) &&
									(inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) &&
									(t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) &&
									(t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) &&
									(t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) &&
										(t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
								{

									//if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible))
									if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) &&
										(b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible))
									{

										ncell_shadow++;
										temp_shadow[inode2] = t.potent[inode1];
										temp_shadow[inode3] = t.potent[inode4];
										temp_shadow[inode6] = t.potent[inode5];
										temp_shadow[inode7] = t.potent[inode8];
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[j_4][inode2] = t.total_deformation[j_4][inode1];
											total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode4];
											total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode5];
											total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode8];
										}
									}
								}
								else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) &&
									(inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) &&
									(t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) &&
									(t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) &&
									(t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) &&
										(t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
								{

									//if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) 
									if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) &&
										(b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible))
									{

										ncell_shadow++;
										temp_shadow[inode4] = t.potent[inode1];
										temp_shadow[inode3] = t.potent[inode2];
										temp_shadow[inode8] = t.potent[inode5];
										temp_shadow[inode7] = t.potent[inode6];
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[j_4][inode4] = t.total_deformation[j_4][inode1];
											total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode2];
											total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode5];
											total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode6];
										}
									}
								}
							}
							
						}
						
                     }
					else {
						

						ncell_shadow = 0;
						for (i = 0; i < ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = t.database.nvtxcell[0][i] - 1;
							inode2 = t.database.nvtxcell[1][i] - 1;
							inode3 = t.database.nvtxcell[2][i] - 1;
							inode4 = t.database.nvtxcell[3][i] - 1;
							inode5 = t.database.nvtxcell[4][i] - 1;
							inode6 = t.database.nvtxcell[5][i] - 1;
							inode7 = t.database.nvtxcell[6][i] - 1;
							inode8 = t.database.nvtxcell[7][i] - 1;

							/*
							TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
							//TOCHKA pall;
							center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
							center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
							center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
							center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
							center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
							center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
							center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
							center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
							*/
							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
							/*
							in_model_temp(p1, ib1, b, lb);
							in_model_temp(p2, ib2, b, lb);
							in_model_temp(p3, ib3, b, lb);
							in_model_temp(p4, ib4, b, lb);
							in_model_temp(p5, ib5, b, lb);
							in_model_temp(p6, ib6, b, lb);
							in_model_temp(p7, ib7, b, lb);
							in_model_temp(p8, ib8, b, lb);
							*/

							ib1 = t.whot_is_block[inode1];
							ib2 = t.whot_is_block[inode2];
							ib3 = t.whot_is_block[inode3];
							ib4 = t.whot_is_block[inode4];
							ib5 = t.whot_is_block[inode5];
							ib6 = t.whot_is_block[inode6];
							ib7 = t.whot_is_block[inode7];
							ib8 = t.whot_is_block[inode8];


							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) &&
								(b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
								ncell_shadow++;
							}
						}
						//ncell_shadow = ncell;

						//printf("debug %lld\n", ncell_shadow); getchar(); // debug
					}
					
                }
				if (ncell_shadow == 0) {
					// У нас совсем нету твёрдого тела, поэтому мы показываем только  жидкость.
					ncell_shadow = ncell;
					ionly_solid_visible = 0;
				}
				if (1 && 2 != steady_or_unsteady_global_determinant) {
					// Полный набор искомых величин и теплопередача и гидродинамика:
					
						//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					    if (f[0].iflowregime == RANS_MENTER_SST) {
							fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, Re_y, k, omega, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
					    }
					    else {
						    fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, Re_y, k, eps, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
					    }
				}
				else {
					// 19,03,2019
					fprintf(fp, "\nVARIABLES = x, y, z, quolity, log(quolity)\n");// Только вызов Мешера
				}

#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif




				if (bvery_big_memory) {
					if (lite_export) {
						// запись x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						if (1 && (2 == steady_or_unsteady_global_determinant)) {
							// quolity
							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
								volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.6f ", dmax/dmin);
								if (i % 10 == 0) fprintf(fp, "\n");
							}

							// log10(quolity)
							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
								volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.6f ", log10(dmax / dmin));
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
					else {
						// запись x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						if (1 && (2 == steady_or_unsteady_global_determinant)) {
							// quolity
							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
								volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.16f ", dmax / dmin);
								if (i % 10 == 0) fprintf(fp, "\n");
							}

							// log10(quolity)
							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
								volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.16f ", log10(dmax / dmin));
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части

		// Запись поля температур производится всегда.
		if (1 && 2 != steady_or_unsteady_global_determinant) {

		// запись температуры
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == 1) {
					fprintf(fp, "%+.3f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.3f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.3f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.3f ", t.potent[i]);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}
		else {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == 1) {
					fprintf(fp, "%+.16f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.16f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.16f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.16f ", t.potent[i]);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}




		fprintf(fp, "\n");

		// Lam
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.4f ", t.prop[LAM][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < t.maxbound; i++) {
					fprintf(fp, "%+.4f ", t.prop_b[LAM][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}
		else {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", t.prop[LAM][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < t.maxbound; i++) {
					fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		fprintf(fp, "\n");

		// Запись гидродинамических величин если необходимо:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx*svx + svy*svy + svz*svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VX][i + maxelm] * f[idfluid].potent[VX][i + maxelm];
					doublereal svy = f[idfluid].potent[VY][i + maxelm] * f[idfluid].potent[VY][i + maxelm];
					doublereal svz = f[idfluid].potent[VZ][i + maxelm] * f[idfluid].potent[VZ][i + maxelm];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}


			fprintf(fp, "\n");

			// Pressure
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Pressure
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// PAM
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// PAM
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				}
				else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Rho
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mu
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // Турбулентная динамическая вязкость
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);

					//fprintf(fp, "%+.16f ", 1.0*f[t.ptr[1][i]].icolor_different_fluid_domain[t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// MUT
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					//fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm] / f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST)||
						(f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime == ZEROEQMOD) || 
						(f[0].iflowregime == SMAGORINSKY)||
						(f[0].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == RANS_MENTER_SST) ||
						(f[0].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // Завихрённость - модуль ротора скорости.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Curl
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i + maxelm]); // GRADZVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Re_y
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					if ((f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)||(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST)) {
						doublereal speed_or_sqrt_k = 0.0;
						if (f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][t.ptr[0][i]]));
						}
						if (f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY][t.ptr[0][i]]));
						}
						//speed_or_sqrt_k = sqrt(f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] + f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] + f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
						integer id_loc = t.ptr[0][i];
						doublereal Re_y = 0.0;
						if (id_loc < f[t.ptr[1][i]].maxelm) {
							Re_y = speed_or_sqrt_k *
								f[t.ptr[1][i]].rdistWall[t.ptr[0][i]] * f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU][t.ptr[0][i]];
						}
						else {
							Re_y = speed_or_sqrt_k *
								f[t.ptr[1][i]].rdistWall[t.ptr[0][i]] * f[t.ptr[1][i]].prop_b[RHO][t.ptr[0][i] - f[t.ptr[1][i]].maxelm] / f[t.ptr[1][i]].prop_b[MU][t.ptr[0][i] - f[t.ptr[1][i]].maxelm];

						}
						fprintf(fp, "%+.16f ", Re_y);
					}
					else {
						fprintf(fp, "%+.16f ", 0.0);
					}
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Re_y
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if ((f[idfluid].iflowregime == RANS_STANDART_K_EPS)||(f[idfluid].iflowregime == RANS_MENTER_SST)) {
						doublereal speed_or_sqrt_k = 0.0;
						if (f[idfluid].iflowregime == RANS_STANDART_K_EPS) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[idfluid].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i + maxelm]));
						}
						if (f[idfluid].iflowregime == RANS_MENTER_SST) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[idfluid].potent[TURBULENT_KINETIK_ENERGY][i + maxelm]));
						}
						//speed_or_sqrt_k = sqrt(f[idfluid].potent[VX][i+ maxelm] * f[idfluid].potent[VX][i+ maxelm] + f[idfluid].potent[VY][i+ maxelm] * f[idfluid].potent[VY][i+ maxelm] + f[idfluid].potent[VZ][i+ maxelm] * f[idfluid].potent[VZ][i+ maxelm]);
						doublereal Re_y = speed_or_sqrt_k *
							f[idfluid].rdistWall[i + maxelm] * f[idfluid].prop_b[RHO][i] / f[idfluid].prop_b[MU][i];

						fprintf(fp, "%+.16f ", Re_y); // Re_y
					}
					else {
						fprintf(fp, "%+.16f ", 0.0);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// TURBULENT_KINETIK_ENERGY_STD_K_EPS
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					if (f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][t.ptr[0][i]]);
					}
					else if (f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// TURBULENT_KINETIK_ENERGY_STD_K_EPS
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if (f[idfluid].iflowregime == RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i + maxelm]); // TURBULENT_KINETIK_ENERGY_STD_K_EPS
					}
					else if (f[idfluid].iflowregime == RANS_MENTER_SST) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[TURBULENT_KINETIK_ENERGY][i + maxelm]); // TURBULENT_KINETIK_ENERGY
					}
					else {
						fprintf(fp, "%+.16f ", 0.0);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					if (f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][t.ptr[0][i]]);
					}
					else if (f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][t.ptr[0][i]]);
					}
					else {
						fprintf(fp, "%+.16f ", 0.0);
					}
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if (f[idfluid].iflowregime == RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i + maxelm]); // TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS
					}
					else if (f[idfluid].iflowregime == RANS_MENTER_SST) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i + maxelm]); // TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA
					}
					else {
						fprintf(fp, "%+.16f ", 0.0);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


		}

		doublereal *Tx = nullptr;
		doublereal *Ty = nullptr;
		doublereal *Tz = nullptr;
		Tx = new doublereal[t.maxelm + t.maxbound];
		Ty = new doublereal[t.maxelm + t.maxbound];
		Tz = new doublereal[t.maxelm + t.maxbound];

		// инициализация нулём.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// нахождение градиентов.
		for (i = 0; i<t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, false,
				t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
		}

		for (i = 0; i<t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, true,
				t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
		}


		// Сохранение в файл.

		if (lite_export) {
			doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

			// Heat Flux X
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Tx[i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Tx[i]);

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tx[i]);
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Heat Flux Y
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Ty[i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Ty[i]);

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Heat Flux Z
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Tz[i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Tz[i]);

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


			// Mag Heat Flux
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));


					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {
						i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


			// log10 Mag Heat Flux
			const doublereal eps_min = 2.0;
			const doublereal d_stub = 0.0;
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));

					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					if (buf1 > eps_min) {
						buf1 = log10(buf1);
					}
					else {
						buf1 = d_stub;
					}
					if (buf2 > eps_min) {
						buf2 = log10(buf2);
					}
					else {
						buf2 = d_stub;
					}
					if (buf3 > eps_min) {
						buf3 = log10(buf3);
					}
					else {
						buf3 = d_stub;
					}
					if (buf4 > eps_min) {
						buf4 = log10(buf4);
					}
					else {
						buf4 = d_stub;
					}
					if (buf5 > eps_min) {
						buf5 = log10(buf5);
					}
					else {
						buf5 = d_stub;
					}
					if (buf6 > eps_min) {
						buf6 = log10(buf6);
					}
					else {
						buf6 = d_stub;
					}
					if (buf7 > eps_min) {
						buf7 = log10(buf7);
					}
					else {
						buf7 = d_stub;
					}
					if (buf8 > eps_min) {
						buf8 = log10(buf8);
					}
					else {
						buf8 = d_stub;
					}
					if (buf9 > eps_min) {
						buf9 = log10(buf9);
					}
					else {
						buf9 = d_stub;
					}

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {
						i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}

						//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

		}
		else {
			doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

			// Heat Flux X
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 1;
					buf1 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 2;
					buf2 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 3;
					buf3 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 4;
					buf4 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 5;
					buf5 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 6;
					buf6 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 7;
					buf7 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 8;
					buf8 = -t.prop[LAM][i] * Tx[i];
					i = i1 + 9;
					buf9 = -t.prop[LAM][i] * Tx[i];

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = -t.prop[LAM][i] * Tx[i];
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 1;
						buf1 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 2;
						buf2 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 3;
						buf3 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 4;
						buf4 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 5;
						buf5 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 6;
						buf6 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 7;
						buf7 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 8;
						buf8 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						i = i1 + 9;
						buf9 = -t.prop_b[LAM][i] * Tx[i + maxelm];


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop_b[LAM][i] * Tx[i + maxelm];
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Heat Flux Y
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 1;
					buf1 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 2;
					buf2 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 3;
					buf3 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 4;
					buf4 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 5;
					buf5 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 6;
					buf6 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 7;
					buf7 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 8;
					buf8 = -t.prop[LAM][i] * Ty[i];
					i = i1 + 9;
					buf9 = -t.prop[LAM][i] * Ty[i];

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = -t.prop[LAM][i] * Ty[i];
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 1;
						buf1 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 2;
						buf2 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 3;
						buf3 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 4;
						buf4 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 5;
						buf5 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 6;
						buf6 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 7;
						buf7 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 8;
						buf8 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						i = i1 + 9;
						buf9 = -t.prop_b[LAM][i] * Ty[i + maxelm];


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop_b[LAM][i] * Ty[i + maxelm];
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Heat Flux Z
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 1;
					buf1 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 2;
					buf2 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 3;
					buf3 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 4;
					buf4 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 5;
					buf5 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 6;
					buf6 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 7;
					buf7 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 8;
					buf8 = -t.prop[LAM][i] * Tz[i];
					i = i1 + 9;
					buf9 = -t.prop[LAM][i] * Tz[i];

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = -t.prop[LAM][i] * Tz[i];
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 1;
						buf1 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 2;
						buf2 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 3;
						buf3 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 4;
						buf4 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 5;
						buf5 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 6;
						buf6 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 7;
						buf7 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 8;
						buf8 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						i = i1 + 9;
						buf9 = -t.prop_b[LAM][i] * Tz[i + maxelm];


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop_b[LAM][i] * Tz[i + maxelm];
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


			// Mag Heat Flux
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));


					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {
						i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.16f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


			// log10 Mag Heat Flux
			const doublereal eps_min = 2.0;
			const doublereal d_stub = 0.0;
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));

					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					if (buf1 > eps_min) {
						buf1 = log10(buf1);
					}
					else {
						buf1 = d_stub;
					}
					if (buf2 > eps_min) {
						buf2 = log10(buf2);
					}
					else {
						buf2 = d_stub;
					}
					if (buf3 > eps_min) {
						buf3 = log10(buf3);
					}
					else {
						buf3 = d_stub;
					}
					if (buf4 > eps_min) {
						buf4 = log10(buf4);
					}
					else {
						buf4 = d_stub;
					}
					if (buf5 > eps_min) {
						buf5 = log10(buf5);
					}
					else {
						buf5 = d_stub;
					}
					if (buf6 > eps_min) {
						buf6 = log10(buf6);
					}
					else {
						buf6 = d_stub;
					}
					if (buf7 > eps_min) {
						buf7 = log10(buf7);
					}
					else {
						buf7 = d_stub;
					}
					if (buf8 > eps_min) {
						buf8 = log10(buf8);
					}
					else {
						buf8 = d_stub;
					}
					if (buf9 > eps_min) {
						buf9 = log10(buf9);
					}
					else {
						buf9 = d_stub;
					}

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {
						i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}

						//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.16f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


		}


		// Освобождение оперативной памяти.
		if (Tx != nullptr) {
			delete[] Tx;
		}
		if (Ty != nullptr) {
			delete[] Ty;
		}
		if (Tz != nullptr) {
			delete[] Tz;
		}

		fprintf(fp, "\n");

		for (integer j_6 = 0; j_6 < 4; j_6++) {

			// запись полной деформации
			if (lite_export) {
				for (i = 0; i < maxelm; i++) {


					if (ionly_solid_visible == 1) {
						fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
						//printf("%e \n", total_deformation_shadow[j_6][i]);
					}
					else {
						fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						//printf("%e \n", t.total_deformation[j_6][i]);
					}
					if (i % 10 == 0) {
						fprintf(fp, "\n");
						//getchar();
					}
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
			else {
				for (i = 0; i < maxelm; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
					}
					else {
						fprintf(fp, "%e ", t.total_deformation[j_6][i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}


			fprintf(fp, "\n");

		}
		} // if (1 && 2 != steady_or_unsteady_global_determinant) 

		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					/*
					TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
					*/
					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					/*
					in_model_temp(p1, ib1, b, lb);
					in_model_temp(p2, ib2, b, lb);
					in_model_temp(p3, ib3, b, lb);
					in_model_temp(p4, ib4, b, lb);
					in_model_temp(p5, ib5, b, lb);
					in_model_temp(p6, ib6, b, lb);
					in_model_temp(p7, ib7, b, lb);
					in_model_temp(p8, ib8, b, lb);
					*/

					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];


					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
						integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
						integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
						integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


						integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
						integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
						integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
						integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



						integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
						integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
						integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
						integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;


						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
						pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
						pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);

						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
					}
					else {
					
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);

						if ((ib1 != t.whot_is_block[inode1]) || (ib2 != t.whot_is_block[inode2]) || (ib3 != t.whot_is_block[inode3]) || (ib4 != t.whot_is_block[inode4]) ||
							(ib5 != t.whot_is_block[inode5]) || (ib6 != t.whot_is_block[inode6]) || (ib7 != t.whot_is_block[inode7]) || (ib8 != t.whot_is_block[inode8])) {
							printf("SIGNAL MESSAGE ZAMENA NEVERNA."); getchar();
						}
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];

						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}
		}
		else {

		
#ifdef MINGW_COMPILLER
		err = 0;
		fp1=fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
		if (fp1==NULL) err = 1;
#else
		err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	if (temp_shadow != nullptr) {
		delete[] temp_shadow;
		temp_shadow = nullptr;
	}
	//total_deformation

	if (total_deformation_shadow != nullptr) {
		for (integer j_6 = 0; j_6 < 4; j_6++) {
			if (total_deformation_shadow[j_6] != nullptr) {
				delete[] total_deformation_shadow[j_6];
				total_deformation_shadow[j_6] = nullptr;
			}
		}
		delete[] total_deformation_shadow;
		total_deformation_shadow = nullptr;
	}


	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);



} // exporttecplotxy360T_3D_part2

void export_tecplot_temperature_ass(integer** &nvtx, bool* &bcheck_visible, TOCHKA* &pa, doublereal* &potent,
	doublereal* &lam_for_export, doublereal* &Txgl, doublereal* &Tygl, doublereal* &Tzgl, doublereal* &HeatFluxMaggl,
	 integer maxelm, integer ncell) {
	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	err = 0;
	fp=fopen64("ALICEFLOW0_08_temp.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "ALICEFLOW0_08_temp.PLT", "wb");
#endif
	

	if ((err) != 0) {
		printf("Create File temp Error in function export_tecplot_temperature_ass in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_08\"\n");

		// запись имён переменных
		fprintf(fp, "VARIABLES = x, y, z, T, lam, Tx, Ty, Tz, magHeatFlux, log10MagHeatFlux \n");  // печатается только поле температур

												   // запись информации о зонах
		integer ncell_actual = 0;
		for (integer i = 0; i < ncell; i++) {
			if (bcheck_visible[i]) {
				ncell_actual++;
			}
		}
		fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_actual);


		// запись x
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", pa[i].x);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись y
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", pa[i].y);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись z
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", pa[i].z);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись Температуры
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", potent[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись Теплопроводности
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", lam_for_export[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		
		// запись -lambda*gradTx
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", Txgl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись -lambda*gradTy
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", Tygl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись -lambda*gradTz
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", Tzgl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись heatfluxMag
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", HeatFluxMaggl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись log10_heatfluxMag
		for (integer i = 0; i < maxelm; i++) {
			if (HeatFluxMaggl[i] > 1.0e-10) {
				fprintf(fp, "%+.6f ", log10(HeatFluxMaggl[i]));
			}
			else {
				fprintf(fp, "%+.6f ", -10.0);
			}
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		for (integer i = 0; i < ncell; i++) {
			if (bcheck_visible[i]) {
				// где то выше по коду перепутанпорядок в nvtx. Здесь он восстановлен.
				if (b_on_adaptive_local_refinement_mesh) {
					fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
				}
				else {
					fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld \n", nvtx[0][i], nvtx[1][i], nvtx[3][i], nvtx[2][i], nvtx[4][i], nvtx[5][i], nvtx[7][i], nvtx[6][i]);
				}
			}
		}

		fclose(fp);
		printf("file succsefull writing\n");

	}
}


void exporttecplot_assembles_mesh(TEMPER &t, integer lu, UNION* &my_union) {
	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");
#endif
	

	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplot_assembles_mesh in my_export_tecplot3.c\n");
		system("pause");

	}
	else {
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

		// запись имён переменных
		fprintf(fp, "VARIABLES = x, y, z \n");  // печатается только поле температур

		

		integer maxelm_global = t.maxnod;
		integer ncell_shadow_gl = t.maxelm;
		for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
			maxelm_global += my_union[iu_74].t.maxnod;			
			ncell_shadow_gl += my_union[iu_74].t.maxelm;
		}


		// запись информации о зонах
		
			fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global , ncell_shadow_gl);
		
			TOCHKA* pa_gl = new TOCHKA[maxelm_global];
			integer pa_gl_count = -1;
			/*
			for (integer i = 0; i < t.maxnod; i++) {
				pa_gl_count++;
				pa_gl[pa_gl_count] = t.pa[i];
			}
			const doublereal eps = 1.0e-23;
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					bool bfound = false;
					for (integer j = 0; j <= pa_gl_count; j++) {
						if ((fabs(pa_gl[j].x - my_union[iunion_scan].t.pa[i].x) < eps)&&
							(fabs(pa_gl[j].y - my_union[iunion_scan].t.pa[i].y) < eps)&&
							(fabs(pa_gl[j].z - my_union[iunion_scan].t.pa[i].z) < eps)) {
							bfound = true;
							break;
						}
					}
					if (!bfound) {
						pa_gl_count++;
						pa_gl[pa_gl_count] = my_union[iunion_scan].t.pa[i];
					}
				}
			}
			printf("maxnod=%d maxnod_path=%d compare=%d\n", maxelm_global, pa_gl_count+1,abs(maxelm_global-(pa_gl_count + 1)));
			*/

			// запись x
			for (integer i = 0; i < t.maxnod; i++) {
				fprintf(fp, "%+.6f ", t.pa[i].x);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// запись x
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					fprintf(fp, "%+.6f ", my_union[iunion_scan].t.pa[i].x);
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}
			

			fprintf(fp, "\n");

			// запись y
			for (integer i = 0; i < t.maxnod; i++) {
				fprintf(fp, "%+.6f ", t.pa[i].y);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// запись x
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					fprintf(fp, "%+.6f ", my_union[iunion_scan].t.pa[i].y);
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}
			

			fprintf(fp, "\n");

			// запись z
			for (integer i = 0; i < t.maxnod; i++) {
				fprintf(fp, "%+.6f ", t.pa[i].z);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// запись z
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					fprintf(fp, "%+.6f ", my_union[iunion_scan].t.pa[i].z);
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}
			
			fprintf(fp, "\n");

			integer iadd = 0;

			for (integer i = 0; i < t.maxelm; i++) {
				fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld \n", t.nvtx[0][i], t.nvtx[1][i], t.nvtx[3][i], t.nvtx[2][i], t.nvtx[4][i], t.nvtx[5][i], t.nvtx[7][i], t.nvtx[6][i]);
			}
			
			iadd += t.maxnod;

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (integer i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld \n", iadd+my_union[iunion_scan].t.nvtx[0][i], iadd + my_union[iunion_scan].t.nvtx[1][i], iadd + my_union[iunion_scan].t.nvtx[3][i], iadd + my_union[iunion_scan].t.nvtx[2][i], iadd + my_union[iunion_scan].t.nvtx[4][i], iadd + my_union[iunion_scan].t.nvtx[5][i], iadd + my_union[iunion_scan].t.nvtx[7][i], iadd + my_union[iunion_scan].t.nvtx[6][i]);
				}
				iadd += my_union[iunion_scan].t.maxnod;
			}

			if (pa_gl != nullptr) {
				delete[] pa_gl;
			}

			fclose(fp);

	}

}

  // 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
  // а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
  // проверка построеной сетки
  // экспорт результата расчёта в программу tecplot360
  // часть 2.
void exporttecplotxy360T_3D_part2_assembles(integer maxelm, integer ncell, 
	FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, 
	bool bextendedprint, integer ikey, integer lu, UNION* &my_union,
	BLOCK* &b, integer &lb)
{
	const bool lite_export = true;
	// 16 знаков после запятой сохранять никому ненужно,
	// вполне достаточно шести знаков.

#if doubleintprecision == 1
	printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
#else
	printf("ionly_solid_visible =%d\n", ionly_solid_visible);
#endif


	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal** temp_shadow = nullptr;
	temp_shadow = new doublereal*[lu+1];
	if (ionly_solid_visible == 1) {
		temp_shadow[0] = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[0][i_1] = t.potent[i_1];
		}
	}
	for (integer i = 0; i < lu; i++) {
		if (ionly_solid_visible == 1) {
			temp_shadow[i+1] = new doublereal[my_union[i].t.maxelm + my_union[i].t.maxbound];
			for (integer i_1 = 0; i_1 < my_union[i].t.maxelm + my_union[i].t.maxbound; i_1++) {
				temp_shadow[i+1][i_1] = my_union[i].t.potent[i_1];
			}
		}
	}

	doublereal*** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == 1) {
	total_deformation_shadow = new doublereal**[lu + 1];
		total_deformation_shadow[0] = new doublereal*[4];
		for (integer j_1 = 0; j_1 < 4; j_1++) {
			total_deformation_shadow[0][j_1] = new doublereal[t.maxelm + t.maxbound];
			for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
				total_deformation_shadow[0][j_1][i_1] = t.total_deformation[j_1][i_1];
			}
		}
	}
	for (integer i = 0; i < lu; i++) {
		if (ionly_solid_visible == 1) {
			total_deformation_shadow[i+1] = new doublereal*[4];
			for (integer j_1 = 0; j_1 < 4; j_1++) {
				total_deformation_shadow[i+1][j_1] = new doublereal[my_union[i].t.maxelm + my_union[i].t.maxbound];
				for (integer i_1 = 0; i_1 < my_union[i].t.maxelm + my_union[i].t.maxbound; i_1++) {
					total_deformation_shadow[i+1][j_1][i_1] = my_union[i].t.total_deformation[j_1][i_1];
				}
			}
		}
	}

	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	err = 0;
	switch (ikey) {
	case 0: fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1: fp = fopen64("ALICEFLOW0_27_temp.PLT", "wb");  break; // То что нужно для отчёта Алексею.
	default: fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
}
	if (fp==NULL) err = 1;
#else
	switch (ikey) {
	case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp.PLT", "wb");  break; // То что нужно для отчёта Алексею.
	default: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	}
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2_assembles in file my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1=fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1==NULL) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		if (bOk)
		{


			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // считаем что температура всегда активна
				}
			}

			if (ivarexport == 1) {

				integer* ncell_shadow = nullptr;
				ncell_shadow = new integer[lu + 1];
				for (i = 0; i < lu + 1; i++) ncell_shadow[i] = 0;


				for (i = 0; i < ncell; i++) {

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					// Использование быстродействующей хеш таблицы whot_is_block значительно
					// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
					ib1 = t.whot_is_block[inode1];
					//in_model_temp(p1, ib1, b, lb);	
					//if (ib1 != t.whot_is_block[inode1]) {
					//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
					//getchar();
					//}
					ib2 = t.whot_is_block[inode2];
					//in_model_temp(p2, ib2, b, lb);
					//if (ib2 != t.whot_is_block[inode2]) {
					//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
					//getchar();
					//}
					ib3 = t.whot_is_block[inode3];
					//in_model_temp(p3, ib3, b, lb);
					//if (ib3 != t.whot_is_block[inode3]) {
					//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
					//getchar();
					//}
					ib4 = t.whot_is_block[inode4];
					//in_model_temp(p4, ib4, b, lb);
					//if (ib4 != t.whot_is_block[inode4]) {
					//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
					//getchar();
					//}
					ib5 = t.whot_is_block[inode5];
					//in_model_temp(p5, ib5, b, lb);
					//if (ib5 != t.whot_is_block[inode5]) {
					//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
					//getchar();
					//}
					ib6 = t.whot_is_block[inode6];
					//in_model_temp(p6, ib6, b, lb);
					//if (ib6 != t.whot_is_block[inode6]) {
					//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
					//getchar();
					//}
					ib7 = t.whot_is_block[inode7];
					//in_model_temp(p7, ib7, b, lb);
					//if (ib7 != t.whot_is_block[inode7]) {
					//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
					//getchar();
					//}
					ib8 = t.whot_is_block[inode8];
					//in_model_temp(p8, ib8, b, lb);
					//if (ib8 != t.whot_is_block[inode8]) {
					//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
					//getchar();
					//}


					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
						ncell_shadow[0]++;
					}
				}

				for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
					for (i = 0; i < my_union[iunion_scan].t.ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = my_union[iunion_scan].t.database.nvtxcell[0][i] - 1;
						inode2 = my_union[iunion_scan].t.database.nvtxcell[1][i] - 1;
						inode3 = my_union[iunion_scan].t.database.nvtxcell[2][i] - 1;
						inode4 = my_union[iunion_scan].t.database.nvtxcell[3][i] - 1;
						inode5 = my_union[iunion_scan].t.database.nvtxcell[4][i] - 1;
						inode6 = my_union[iunion_scan].t.database.nvtxcell[5][i] - 1;
						inode7 = my_union[iunion_scan].t.database.nvtxcell[6][i] - 1;
						inode8 = my_union[iunion_scan].t.database.nvtxcell[7][i] - 1;

						//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						// Использование быстродействующей хеш таблицы whot_is_block значительно
						// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
						ib1 = my_union[iunion_scan].t.whot_is_block[inode1];
						//in_model_temp(p1, ib1, b, lb);	
						//if (ib1 != t.whot_is_block[inode1]) {
						//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
						//getchar();
						//}
						ib2 = my_union[iunion_scan].t.whot_is_block[inode2];
						//in_model_temp(p2, ib2, b, lb);
						//if (ib2 != t.whot_is_block[inode2]) {
						//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
						//getchar();
						//}
						ib3 = my_union[iunion_scan].t.whot_is_block[inode3];
						//in_model_temp(p3, ib3, b, lb);
						//if (ib3 != t.whot_is_block[inode3]) {
						//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
						//getchar();
						//}
						ib4 = my_union[iunion_scan].t.whot_is_block[inode4];
						//in_model_temp(p4, ib4, b, lb);
						//if (ib4 != t.whot_is_block[inode4]) {
						//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
						//getchar();
						//}
						ib5 = my_union[iunion_scan].t.whot_is_block[inode5];
						//in_model_temp(p5, ib5, b, lb);
						//if (ib5 != t.whot_is_block[inode5]) {
						//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
						//getchar();
						//}
						ib6 = my_union[iunion_scan].t.whot_is_block[inode6];
						//in_model_temp(p6, ib6, b, lb);
						//if (ib6 != t.whot_is_block[inode6]) {
						//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
						//getchar();
						//}
						ib7 = my_union[iunion_scan].t.whot_is_block[inode7];
						//in_model_temp(p7, ib7, b, lb);
						//if (ib7 != t.whot_is_block[inode7]) {
						//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
						//getchar();
						//}
						ib8 = my_union[iunion_scan].t.whot_is_block[inode8];
						//in_model_temp(p8, ib8, b, lb);
						//if (ib8 != t.whot_is_block[inode8]) {
						//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
						//getchar();
						//}


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
							ncell_shadow[iunion_scan+1]++;
						}
					}
				}
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// запись имён переменных
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
				fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // печатается только поле температур

				integer maxelm_global = maxelm;
				integer maxbound_global = t.maxbound;
				integer ncell_shadow_gl = ncell_shadow[0];
				for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
					maxelm_global += my_union[iu_74].t.maxelm;
					maxbound_global += my_union[iu_74].t.maxbound;
					ncell_shadow_gl += ncell_shadow[iu_74+1];
				}

#if doubleintprecision == 1
																																																							   // запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#else
																																																							   // запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#endif



				if (bvery_big_memory) {
					// extended print не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
						// запись x
						for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
						// запись y
						for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
						// запись z
						for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			
				delete[] ncell_shadow;
            }
			else if (ivarexport == 3) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer* ncell_shadow = nullptr;
				ncell_shadow = new integer[lu + 1];
				for (i = 0; i < lu + 1; i++) ncell_shadow[i] = 0;

				if (bsolid_static_only) {


					for (i = 0; i < ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						/*
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);
						*/

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
							ncell_shadow[0]++;
						}
					}

					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {

						for (i = 0; i < my_union[iunion_scan].t.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = my_union[iunion_scan].t.database.nvtxcell[0][i] - 1;
							inode2 = my_union[iunion_scan].t.database.nvtxcell[1][i] - 1;
							inode3 = my_union[iunion_scan].t.database.nvtxcell[2][i] - 1;
							inode4 = my_union[iunion_scan].t.database.nvtxcell[3][i] - 1;
							inode5 = my_union[iunion_scan].t.database.nvtxcell[4][i] - 1;
							inode6 = my_union[iunion_scan].t.database.nvtxcell[5][i] - 1;
							inode7 = my_union[iunion_scan].t.database.nvtxcell[6][i] - 1;
							inode8 = my_union[iunion_scan].t.database.nvtxcell[7][i] - 1;

							/*
							TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
							//TOCHKA pall;
							center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
							center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
							center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
							center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
							center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
							center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
							center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
							center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
							*/
							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
							/*
							in_model_temp(p1, ib1, b, lb);
							in_model_temp(p2, ib2, b, lb);
							in_model_temp(p3, ib3, b, lb);
							in_model_temp(p4, ib4, b, lb);
							in_model_temp(p5, ib5, b, lb);
							in_model_temp(p6, ib6, b, lb);
							in_model_temp(p7, ib7, b, lb);
							in_model_temp(p8, ib8, b, lb);
							*/

							ib1 = my_union[iunion_scan].t.whot_is_block[inode1];
							ib2 = my_union[iunion_scan].t.whot_is_block[inode2];
							ib3 = my_union[iunion_scan].t.whot_is_block[inode3];
							ib4 = my_union[iunion_scan].t.whot_is_block[inode4];
							ib5 = my_union[iunion_scan].t.whot_is_block[inode5];
							ib6 = my_union[iunion_scan].t.whot_is_block[inode6];
							ib7 = my_union[iunion_scan].t.whot_is_block[inode7];
							ib8 = my_union[iunion_scan].t.whot_is_block[inode8];


							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
								ncell_shadow[iunion_scan+1]++;
							}
						}
					}

					//ncell_shadow = ncell;
				}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{

						for (i = 0; i < t.database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = t.database.nvtxcell[0][i] - 1;
							inode2 = t.database.nvtxcell[1][i] - 1;
							inode3 = t.database.nvtxcell[2][i] - 1;
							inode4 = t.database.nvtxcell[3][i] - 1;
							inode5 = t.database.nvtxcell[4][i] - 1;
							inode6 = t.database.nvtxcell[5][i] - 1;
							inode7 = t.database.nvtxcell[6][i] - 1;
							inode8 = t.database.nvtxcell[7][i] - 1;

							integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
							integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
							integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
							integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


							integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
							integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
							integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
							integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



							integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
							integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
							integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
							integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;

							TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
							center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
							center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
							center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
							center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
							center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
							center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
							center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
							center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
							pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
							pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
							pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
							in_model_temp(p1, ib1, b, lb);
							in_model_temp(p2, ib2, b, lb);
							in_model_temp(p3, ib3, b, lb);
							in_model_temp(p4, ib4, b, lb);
							in_model_temp(p5, ib5, b, lb);
							in_model_temp(p6, ib6, b, lb);
							in_model_temp(p7, ib7, b, lb);
							in_model_temp(p8, ib8, b, lb);

							// визуализация только твёрдого тела.
							// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
							if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
								&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
							{
								if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
									ncell_shadow[0]++;
								}

							}
							else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
							{

								if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
									ncell_shadow[0]++;
									temp_shadow[0][inode5] = t.potent[inode1];
									temp_shadow[0][inode6] = t.potent[inode2];
									temp_shadow[0][inode7] = t.potent[inode3];
									temp_shadow[0][inode8] = t.potent[inode4];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[0][j_4][inode5] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[0][j_4][inode6] = t.total_deformation[j_4][inode2];
										total_deformation_shadow[0][j_4][inode7] = t.total_deformation[j_4][inode3];
										total_deformation_shadow[0][j_4][inode8] = t.total_deformation[j_4][inode4];
									}
								}
							}
							else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
							{

								if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

									ncell_shadow[0]++;
									temp_shadow[0][inode2] = t.potent[inode1];
									temp_shadow[0][inode3] = t.potent[inode4];
									temp_shadow[0][inode6] = t.potent[inode5];
									temp_shadow[0][inode7] = t.potent[inode8];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[0][j_4][inode2] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[0][j_4][inode3] = t.total_deformation[j_4][inode4];
										total_deformation_shadow[0][j_4][inode6] = t.total_deformation[j_4][inode5];
										total_deformation_shadow[0][j_4][inode7] = t.total_deformation[j_4][inode8];
									}
								}
							}
							else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
							{

								if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

									ncell_shadow[0]++;
									temp_shadow[0][inode4] = t.potent[inode1];
									temp_shadow[0][inode3] = t.potent[inode2];
									temp_shadow[0][inode8] = t.potent[inode5];
									temp_shadow[0][inode7] = t.potent[inode6];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[0][j_4][inode4] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[0][j_4][inode3] = t.total_deformation[j_4][inode2];
										total_deformation_shadow[0][j_4][inode8] = t.total_deformation[j_4][inode5];
										total_deformation_shadow[0][j_4][inode7] = t.total_deformation[j_4][inode6];
									}
								}
							}
						}

						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							for (i = 0; i < my_union[iunion_scan].t.database.ncell; i++) {

								integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
								inode1 = my_union[iunion_scan].t.database.nvtxcell[0][i] - 1;
								inode2 = my_union[iunion_scan].t.database.nvtxcell[1][i] - 1;
								inode3 = my_union[iunion_scan].t.database.nvtxcell[2][i] - 1;
								inode4 = my_union[iunion_scan].t.database.nvtxcell[3][i] - 1;
								inode5 = my_union[iunion_scan].t.database.nvtxcell[4][i] - 1;
								inode6 = my_union[iunion_scan].t.database.nvtxcell[5][i] - 1;
								inode7 = my_union[iunion_scan].t.database.nvtxcell[6][i] - 1;
								inode8 = my_union[iunion_scan].t.database.nvtxcell[7][i] - 1;

								integer inode2W = my_union[iunion_scan].t.sosedi[WSIDE][inode1].iNODE1;
								integer inode3W = my_union[iunion_scan].t.sosedi[WSIDE][inode4].iNODE1;
								integer inode6W = my_union[iunion_scan].t.sosedi[WSIDE][inode5].iNODE1;
								integer inode7W = my_union[iunion_scan].t.sosedi[WSIDE][inode8].iNODE1;


								integer inode5B = my_union[iunion_scan].t.sosedi[BSIDE][inode1].iNODE1;
								integer inode6B = my_union[iunion_scan].t.sosedi[BSIDE][inode2].iNODE1;
								integer inode7B = my_union[iunion_scan].t.sosedi[BSIDE][inode3].iNODE1;
								integer inode8B = my_union[iunion_scan].t.sosedi[BSIDE][inode4].iNODE1;



								integer inode3S = my_union[iunion_scan].t.sosedi[SSIDE][inode2].iNODE1;
								integer inode4S = my_union[iunion_scan].t.sosedi[SSIDE][inode1].iNODE1;
								integer inode7S = my_union[iunion_scan].t.sosedi[SSIDE][inode6].iNODE1;
								integer inode8S = my_union[iunion_scan].t.sosedi[SSIDE][inode5].iNODE1;

								TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
								center_cord3D(inode1, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p1, 100);
								center_cord3D(inode2, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p2, 100);
								center_cord3D(inode3, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p3, 100);
								center_cord3D(inode4, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p4, 100);
								center_cord3D(inode5, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p5, 100);
								center_cord3D(inode6, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p6, 100);
								center_cord3D(inode7, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p7, 100);
								center_cord3D(inode8, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p8, 100);
								pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
								pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
								pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

								integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
								in_model_temp(p1, ib1, b, lb);
								in_model_temp(p2, ib2, b, lb);
								in_model_temp(p3, ib3, b, lb);
								in_model_temp(p4, ib4, b, lb);
								in_model_temp(p5, ib5, b, lb);
								in_model_temp(p6, ib6, b, lb);
								in_model_temp(p7, ib7, b, lb);
								in_model_temp(p8, ib8, b, lb);

								// визуализация только твёрдого тела.
								// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
								if (((my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1)
									&& (my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))
								{
									if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
										ncell_shadow[iunion_scan+1]++;
									}

								}
								else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) && (!((my_union[iunion_scan].t.ptr[1][inode5B] == -1) && (my_union[iunion_scan].t.ptr[1][inode6B] == -1) && (my_union[iunion_scan].t.ptr[1][inode7B] == -1) && (my_union[iunion_scan].t.ptr[1][inode8B] == -1))) && (!((my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))))
								{

									if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
										ncell_shadow[iunion_scan + 1]++;
										temp_shadow[iunion_scan+1][inode5] = my_union[iunion_scan].t.potent[inode1];
										temp_shadow[iunion_scan + 1][inode6] = my_union[iunion_scan].t.potent[inode2];
										temp_shadow[iunion_scan + 1][inode7] = my_union[iunion_scan].t.potent[inode3];
										temp_shadow[iunion_scan + 1][inode8] = my_union[iunion_scan].t.potent[inode4];
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[iunion_scan + 1][j_4][inode5] = my_union[iunion_scan].t.total_deformation[j_4][inode1];
											total_deformation_shadow[iunion_scan + 1][j_4][inode6] = my_union[iunion_scan].t.total_deformation[j_4][inode2];
											total_deformation_shadow[iunion_scan + 1][j_4][inode7] = my_union[iunion_scan].t.total_deformation[j_4][inode3];
											total_deformation_shadow[iunion_scan + 1][j_4][inode8] = my_union[iunion_scan].t.total_deformation[j_4][inode4];
										}
									}
								}
								else if (((inode2W >= 0) && (inode2W < my_union[iunion_scan].t.maxelm) && (inode3W >= 0) && (inode3W < my_union[iunion_scan].t.maxelm) && (inode6W >= 0) && (inode6W < my_union[iunion_scan].t.maxelm) && (inode7W >= 0) && (inode7W < my_union[iunion_scan].t.maxelm) && (my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) && (my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1) && (!((my_union[iunion_scan].t.ptr[1][inode2W] == -1) && (my_union[iunion_scan].t.ptr[1][inode3W] == -1) && (my_union[iunion_scan].t.ptr[1][inode6W] == -1) && (my_union[iunion_scan].t.ptr[1][inode7W] == -1))) && (!((my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1)))))
								{

									if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

										ncell_shadow[iunion_scan + 1]++;
										temp_shadow[iunion_scan+1][inode2] = my_union[iunion_scan].t.potent[inode1];
										temp_shadow[iunion_scan + 1][inode3] = my_union[iunion_scan].t.potent[inode4];
										temp_shadow[iunion_scan + 1][inode6] = my_union[iunion_scan].t.potent[inode5];
										temp_shadow[iunion_scan + 1][inode7] = my_union[iunion_scan].t.potent[inode8];
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[iunion_scan + 1][j_4][inode2] = my_union[iunion_scan].t.total_deformation[j_4][inode1];
											total_deformation_shadow[iunion_scan + 1][j_4][inode3] = my_union[iunion_scan].t.total_deformation[j_4][inode4];
											total_deformation_shadow[iunion_scan + 1][j_4][inode6] = my_union[iunion_scan].t.total_deformation[j_4][inode5];
											total_deformation_shadow[iunion_scan + 1][j_4][inode7] = my_union[iunion_scan].t.total_deformation[j_4][inode8];
										}
									}
								}
								else if (((inode3S >= 0) && (inode3S < my_union[iunion_scan].t.maxelm) && (inode4S >= 0) && (inode4S < my_union[iunion_scan].t.maxelm) && (inode7S >= 0) && (inode7S < my_union[iunion_scan].t.maxelm) && (inode8S >= 0) && (inode8S < my_union[iunion_scan].t.maxelm) && (my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (!((my_union[iunion_scan].t.ptr[1][inode3S] == -1) && (my_union[iunion_scan].t.ptr[1][inode4S] == -1) && (my_union[iunion_scan].t.ptr[1][inode7S] == -1) && (my_union[iunion_scan].t.ptr[1][inode8S] == -1))) && (!((my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))))
								{

									if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

										ncell_shadow[iunion_scan + 1]++;
										temp_shadow[iunion_scan+1][inode4] = my_union[iunion_scan].t.potent[inode1];
										temp_shadow[iunion_scan + 1][inode3] = my_union[iunion_scan].t.potent[inode2];
										temp_shadow[iunion_scan + 1][inode8] = my_union[iunion_scan].t.potent[inode5];
										temp_shadow[iunion_scan + 1][inode7] = my_union[iunion_scan].t.potent[inode6];
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[iunion_scan + 1][j_4][inode4] = my_union[iunion_scan].t.total_deformation[j_4][inode1];
											total_deformation_shadow[iunion_scan + 1][j_4][inode3] = my_union[iunion_scan].t.total_deformation[j_4][inode2];
											total_deformation_shadow[iunion_scan + 1][j_4][inode8] = my_union[iunion_scan].t.total_deformation[j_4][inode5];
											total_deformation_shadow[iunion_scan + 1][j_4][inode7] = my_union[iunion_scan].t.total_deformation[j_4][inode6];
										}
									}
								}
							}
						}

						doublereal iaccsum = 0;
						for (integer i83 = 0; i83 <= lu; i83++) {
							iaccsum += ncell_shadow[i83];
						}
						if (iaccsum == 0) {
							// У нас совсем нету твёрдого тела, поэтому мы показываем только  жидкость.
							ncell_shadow[0] = ncell;
							for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
								ncell_shadow[iunion_scan+1] += my_union[iunion_scan].t.ncell;
							}
							ionly_solid_visible = 0;
						}
					}
				}
				// Полный набор искомых величин и теплопередача и гидродинамика:
				if (bextendedprint) {
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				}
				else {
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				}

				integer maxelm_global = maxelm;
				integer maxbound_global = t.maxbound;
				integer ncell_shadow_gl = ncell_shadow[0];
				for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
					maxelm_global += my_union[iu_74].t.maxelm;
					maxbound_global += my_union[iu_74].t.maxbound;
					ncell_shadow_gl+= ncell_shadow[iu_74+1];
				}

#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#endif

				if (bvery_big_memory) {
					if (lite_export) {
						// запись x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// запись x
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", my_union[iunion_scan].t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// запись y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// запись y
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", my_union[iunion_scan].t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// запись z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// запись z
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", my_union[iunion_scan].t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
					else {
						// запись x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// запись x
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// запись y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// запись y
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// запись z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// запись z
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
				}
				else {
					// Эта ветка кода теперь никогда не используется,
					// Все хранится в оперативной памяти, т.к. запись и 
					// считывание файла слишком медленные.
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}

				delete[] ncell_shadow;
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части

		// Запись поля температур производится всегда.


		// запись температуры
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == 1) {
					fprintf(fp, "%+.3f ", temp_shadow[0][i]);
				}
				else {
					fprintf(fp, "%+.3f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.3f ", temp_shadow[0][i]);
					}
					else {
						fprintf(fp, "%+.3f ", t.potent[i]);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.3f ", temp_shadow[iunion_scan + 1][i]);
					}
					else {
						fprintf(fp, "%+.3f ", my_union[iunion_scan].t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = my_union[iunion_scan].t.maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%+.3f ", temp_shadow[iunion_scan+1][i]);
						}
						else {
							fprintf(fp, "%+.3f ", my_union[iunion_scan].t.potent[i]);
						}
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

		}
		else {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == 1) {
					fprintf(fp, "%+.16f ", temp_shadow[0][i]);
				}
				else {
					fprintf(fp, "%+.16f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.16f ", temp_shadow[0][i]);
					}
					else {
						fprintf(fp, "%+.16f ", t.potent[i]);
					}
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.16f ", temp_shadow[iunion_scan + 1][i]);
					}
					else {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = my_union[iunion_scan].t.maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%+.16f ", temp_shadow[iunion_scan+1][i]);
						}
						else {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].t.potent[i]);
						}
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

		}




		fprintf(fp, "\n");

		// Lam
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.4f ", t.prop[LAM][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < t.maxbound; i++) {
					fprintf(fp, "%+.4f ", t.prop_b[LAM][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (i = 0; i <  my_union[iunion_scan].t.maxelm; i++) {
					fprintf(fp, "%+.4f ", my_union[iunion_scan].t.prop[LAM][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = 0; i < my_union[iunion_scan].t.maxbound; i++) {
						fprintf(fp, "%+.4f ", my_union[iunion_scan].t.prop_b[LAM][i]);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

		}
		else {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", t.prop[LAM][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < t.maxbound; i++) {
					fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					fprintf(fp, "%+.16f ", my_union[iunion_scan].t.prop[LAM][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = 0; i < my_union[iunion_scan].t.maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].t.prop_b[LAM][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		}
		
		fprintf(fp, "\n");

		// Запись гидродинамических величин если необходимо:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VX][i + maxelm] * f[idfluid].potent[VX][i + maxelm];
					doublereal svy = f[idfluid].potent[VY][i + maxelm] * f[idfluid].potent[VY][i + maxelm];
					doublereal svz = f[idfluid].potent[VZ][i + maxelm] * f[idfluid].potent[VZ][i + maxelm];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Speed
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						doublereal svx = my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VX][my_union[iunion_scan].t.ptr[0][i]] * my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VX][my_union[iunion_scan].t.ptr[0][i]];
						doublereal svy = my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VY][my_union[iunion_scan].t.ptr[0][i]] * my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VY][my_union[iunion_scan].t.ptr[0][i]];
						doublereal svz = my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VZ][my_union[iunion_scan].t.ptr[0][i]] * my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VZ][my_union[iunion_scan].t.ptr[0][i]];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}


				if (bextendedprint) {
					// Speed
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						doublereal svx = my_union[iunion_scan].f[idfluid].potent[VX][i + my_union[iunion_scan].t.maxelm] * my_union[iunion_scan].f[idfluid].potent[VX][i + my_union[iunion_scan].t.maxelm];
						doublereal svy = my_union[iunion_scan].f[idfluid].potent[VY][i + my_union[iunion_scan].t.maxelm] * my_union[iunion_scan].f[idfluid].potent[VY][i + my_union[iunion_scan].t.maxelm];
						doublereal svz = my_union[iunion_scan].f[idfluid].potent[VZ][i + my_union[iunion_scan].t.maxelm] * my_union[iunion_scan].f[idfluid].potent[VZ][i + my_union[iunion_scan].t.maxelm];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Pressure
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Pressure
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {

				// Pressure
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[PRESS][my_union[iunion_scan].t.ptr[0][i]]); // PRESSURE
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Pressure
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[PRESS][i + my_union[iunion_scan].t.maxelm]); // PRESSURE
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// PAM
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// PAM
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {

				// PAM
				for (i = 0; i < maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[PAM][my_union[iunion_scan].t.ptr[0][i]]); // PAM
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// PAM
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[PAM][i + my_union[iunion_scan].t.maxelm]); // PRESSURE
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			
			fprintf(fp, "\n");

			// VX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// VX
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VX][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VX
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[VX][i + my_union[iunion_scan].t.maxelm]); // VX
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// VY
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VY][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VY
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[VY][i + my_union[iunion_scan].t.maxelm]); // VY
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// VZ
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VZ][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VZ
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[VZ][i + my_union[iunion_scan].t.maxelm]); // VZ
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				}
				else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Rho
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Rho
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].prop[RHO][my_union[iunion_scan].t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].diag_coef[VX][i]);
					}
					else fprintf(fp, "%+.16f ", my_union[iunion_scan].t.prop[RHO][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Rho
					for (i = 0; i < my_union[iunion_scan].f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[0].prop_b[RHO][i]);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Mu
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Mu
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].prop[MU][my_union[iunion_scan].t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].slau[VX][i].ap);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Mu
					for (i = 0; i < my_union[iunion_scan].f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[0].prop_b[MU][i]);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Mut // Турбулентная динамическая вязкость
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);

					//fprintf(fp, "%+.16f ", 1.0*f[t.ptr[1][i]].icolor_different_fluid_domain[t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// MUT
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					//fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm] / f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Mut // Турбулентная динамическая вязкость
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[MUT][my_union[iunion_scan].t.ptr[0][i]]);
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[MUT][my_union[iunion_scan].t.ptr[0][i]] / my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].prop[MU][my_union[iunion_scan].t.ptr[0][i]]);

						//fprintf(fp, "%+.16f ", 1.0*my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].icolor_different_fluid_domain[my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// MUT
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[MUT][i + my_union[iunion_scan].t.maxelm]); // MUT
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[MUT][i + my_union[iunion_scan].t.maxelm] / my_union[iunion_scan].f[0].prop_b[MU][i]);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

			}

			fprintf(fp, "\n");

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime == ZEROEQMOD) || 
						(f[0].iflowregime == SMAGORINSKY)||
						(f[0].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == RANS_MENTER_SST) ||
						(f[0].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// для отладки график распределения нумерации контрольных объёмов.
				// или распределение расстояния до стенки Distance_Wall.
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", doublereal(i));
						if ((my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == ZEROEQMOD) ||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == SMAGORINSKY)||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES) ||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == RANS_MENTER_SST) ||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].rdistWall[my_union[iunion_scan].t.ptr[0][i]]);
						}
						else fprintf(fp, "%+.16f ", 0.0);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Distance_Wall.
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						if ((my_union[iunion_scan].f[0].iflowregime == ZEROEQMOD) || 
							(my_union[iunion_scan].f[0].iflowregime == SMAGORINSKY)||
							(my_union[iunion_scan].f[0].iflowregime == RANS_SPALART_ALLMARES)||
							(my_union[iunion_scan].f[0].iflowregime == RANS_MENTER_SST) ||
							(my_union[iunion_scan].f[0].iflowregime == RANS_STANDART_K_EPS)) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].rdistWall[i + my_union[iunion_scan].t.maxelm]); // Distance_Wall
						}
						else fprintf(fp, "%+.16f ", 0.0);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");



			// Curl // Завихрённость - модуль ротора скорости.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Curl
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Curl // Завихрённость - модуль ротора скорости.
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[CURL][my_union[iunion_scan].t.ptr[0][i]]); // CURL FBUF
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Curl
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[CURL][i + my_union[iunion_scan].t.maxelm]); // Curl
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
			
			fprintf(fp, "\n");

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADXVX
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADXVX][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVX
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADXVX][i + my_union[iunion_scan].t.maxelm]); // GRADXVX
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADYVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADYVX
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADYVX][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVX
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADYVX][i + my_union[iunion_scan].t.maxelm]); // GRADYVX
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADZVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADZVX
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].potent[GRADZVX][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVX
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADZVX][i + my_union[iunion_scan].t.maxelm]); // GRADZVX
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADXVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADXVY
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADXVY][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVY
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADXVY][i + my_union[iunion_scan].t.maxelm]); // GRADXVY
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADYVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {

				// GRADYVY
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADYVY][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVY
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADYVY][i + my_union[iunion_scan].t.maxelm]); // GRADYVY
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADZVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADZVY
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADZVY][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVY
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADZVY][i + my_union[iunion_scan].t.maxelm]); // GRADZVY
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADXVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADXVZ
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADXVZ][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVZ
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADXVZ][i + my_union[iunion_scan].t.maxelm]); // GRADXVZ
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADYVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADYVZ
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADYVZ][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVZ
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADYVZ][i + my_union[iunion_scan].t.maxelm]); // GRADYVZ
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// GRADZVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i + maxelm]); // GRADZVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// GRADZVZ
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[GRADZVZ][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVZ
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[GRADZVZ][i + my_union[iunion_scan].t.maxelm]); // GRADZVZ
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

		}

		doublereal **Tx = nullptr;
		doublereal **Ty = nullptr;
		doublereal **Tz = nullptr;
		Tx = new doublereal*[lu + 1];
		Ty = new doublereal*[lu + 1];
		Tz = new doublereal*[lu + 1];

		Tx[0] = new doublereal[t.maxelm + t.maxbound];
		Ty[0] = new doublereal[t.maxelm + t.maxbound];
		Tz[0] = new doublereal[t.maxelm + t.maxbound];

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			Tx[iunion_scan+1] = new doublereal[my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound];
			Ty[iunion_scan + 1] = new doublereal[my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound];
			Tz[iunion_scan + 1] = new doublereal[my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound];
		}

		// инициализация нулём.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[0][i] = 0.0;
			Ty[0][i] = 0.0;
			Tz[0][i] = 0.0;
		}

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// инициализация нулём.
			for (i = 0; i<my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
				Tx[iunion_scan + 1][i] = 0.0;
				Ty[iunion_scan + 1][i] = 0.0;
				Tz[iunion_scan + 1][i] = 0.0;
			}
		}

		// нахождение градиентов.
		for (i = 0; i<t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, false,
				t.sosedb, Tx[0], Ty[0], Tz[0], t.ilevel_alice);
		}

		for (i = 0; i<t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, true,
				t.sosedb, Tx[0], Ty[0], Tz[0], t.ilevel_alice);
		}


		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// нахождение градиентов.
			for (i = 0; i<my_union[iunion_scan].t.maxelm; i++) {
				// Только внутренние узлы.
				green_gaussTemperature(i, my_union[iunion_scan].t.potent, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa,
					my_union[iunion_scan].t.sosedi, my_union[iunion_scan].t.maxelm, false,
					my_union[iunion_scan].t.sosedb, Tx[iunion_scan + 1], Ty[iunion_scan + 1], Tz[iunion_scan + 1], 
					my_union[iunion_scan].t.ilevel_alice);
			}

			for (i = 0; i<my_union[iunion_scan].t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, my_union[iunion_scan].t.potent, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa,
					my_union[iunion_scan].t.sosedi, my_union[iunion_scan].t.maxelm, true,
					my_union[iunion_scan].t.sosedb, Tx[iunion_scan + 1], Ty[iunion_scan + 1], Tz[iunion_scan + 1],
					my_union[iunion_scan].t.ilevel_alice);
			}
		}

		// Сохранение в файл.

		if (lite_export) {


			doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

			// Heat Flux X
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Tx[0][i]);

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}


			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

				// Heat Flux X
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {

						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan+1][i]);
						i = i1 + 1;
						buf1 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 2;
						buf2 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 3;
						buf3 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 4;
						buf4 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 5;
						buf5 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 6;
						buf6 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 7;
						buf7 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 8;
						buf8 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						i = i1 + 9;
						buf9 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
						if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {

							i = i1;
							buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			}


			fprintf(fp, "\n");

			// Heat Flux Y
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Ty[0][i]);

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Heat Flux Y
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {

						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 1;
						buf1 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 2;
						buf2 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 3;
						buf3 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 4;
						buf4 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 5;
						buf5 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 6;
						buf6 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 7;
						buf7 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 8;
						buf8 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						i = i1 + 9;
						buf9 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
						if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {

							i = i1;
							buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			}

			fprintf(fp, "\n");

			// Heat Flux Z
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Tz[0][i]);

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Heat Flux Z
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {

						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 1;
						buf1 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 2;
						buf2 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 3;
						buf3 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 4;
						buf4 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 5;
						buf5 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 6;
						buf6 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 7;
						buf7 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 8;
						buf8 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						i = i1 + 9;
						buf9 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
						if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {

							i = i1;
							buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			}

			fprintf(fp, "\n");


			// Mag Heat Flux
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));


					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {
						i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
						//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.6f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Mag Heat Flux
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {
						i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 1;
						buf1 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 2;
						buf2 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 3;
						buf3 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 4;
						buf4 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 5;
						buf5 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 6;
						buf6 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 7;
						buf7 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 8;
						buf8 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 9;
						buf9 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
						if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {
							i = i1;
							buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i +  my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							//fprintf(fp, "%+.16f ", -my_union[iunion_scan].t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			}


			fprintf(fp, "\n");


			// log10 Mag Heat Flux
			const doublereal eps_min = 2.0;
			const doublereal d_stub = 0.0;
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));

					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					if (buf1 > eps_min) {
						buf1 = log10(buf1);
					}
					else {
						buf1 = d_stub;
					}
					if (buf2 > eps_min) {
						buf2 = log10(buf2);
					}
					else {
						buf2 = d_stub;
					}
					if (buf3 > eps_min) {
						buf3 = log10(buf3);
					}
					else {
						buf3 = d_stub;
					}
					if (buf4 > eps_min) {
						buf4 = log10(buf4);
					}
					else {
						buf4 = d_stub;
					}
					if (buf5 > eps_min) {
						buf5 = log10(buf5);
					}
					else {
						buf5 = d_stub;
					}
					if (buf6 > eps_min) {
						buf6 = log10(buf6);
					}
					else {
						buf6 = d_stub;
					}
					if (buf7 > eps_min) {
						buf7 = log10(buf7);
					}
					else {
						buf7 = d_stub;
					}
					if (buf8 > eps_min) {
						buf8 = log10(buf8);
					}
					else {
						buf8 = d_stub;
					}
					if (buf9 > eps_min) {
						buf9 = log10(buf9);
					}
					else {
						buf9 = d_stub;
					}

					fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					fprintf(fp, "%+.6f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {
						i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i +maxelm])*(-t.prop_b[LAM][i] * Tx[0][i +maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i +maxelm])*(-t.prop_b[LAM][i] * Ty[0][i +maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i +maxelm])*(-t.prop_b[LAM][i] * Tz[0][i +maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}

						//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[0][i +maxelm]);
						fprintf(fp, "%+.6f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// log10 Mag Heat Flux
				const doublereal eps_min = 2.0;
				const doublereal d_stub = 0.0;
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {
						i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 1;
						buf1 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 2;
						buf2 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 3;
						buf3 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 4;
						buf4 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 5;
						buf5 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 6;
						buf6 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 7;
						buf7 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 8;
						buf8 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						i = i1 + 9;
						buf9 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
						if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {
							i = i1;
							buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}
							if (buf1 > eps_min) {
								buf1 = log10(buf1);
							}
							else {
								buf1 = d_stub;
							}
							if (buf2 > eps_min) {
								buf2 = log10(buf2);
							}
							else {
								buf2 = d_stub;
							}
							if (buf3 > eps_min) {
								buf3 = log10(buf3);
							}
							else {
								buf3 = d_stub;
							}
							if (buf4 > eps_min) {
								buf4 = log10(buf4);
							}
							else {
								buf4 = d_stub;
							}
							if (buf5 > eps_min) {
								buf5 = log10(buf5);
							}
							else {
								buf5 = d_stub;
							}
							if (buf6 > eps_min) {
								buf6 = log10(buf6);
							}
							else {
								buf6 = d_stub;
							}
							if (buf7 > eps_min) {
								buf7 = log10(buf7);
							}
							else {
								buf7 = d_stub;
							}
							if (buf8 > eps_min) {
								buf8 = log10(buf8);
							}
							else {
								buf8 = d_stub;
							}
							if (buf9 > eps_min) {
								buf9 = log10(buf9);
							}
							else {
								buf9 = d_stub;
							}

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}

							//fprintf(fp, "%+.16f ", -my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
							fprintf(fp, "%+.6f ", buf0);


							if ((i1 + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			}

			fprintf(fp, "\n");

		}
		else {
		

			doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

			// Heat Flux X
			for (integer i1 = 0; i1 < maxelm; i1++) {
				if ((i1 + 10) < maxelm) {

					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop[LAM][i] * Tx[0][i]);

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop[LAM][i] * Tx[0][i]);
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < t.maxbound; i1++) {
					if ((i1 + 10) < t.maxbound) {

						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop_b[LAM][i] * Tx[0][i + maxelm]);
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}


			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {

			buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

			// Heat Flux X
			for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
				if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {

					i = i1;
					buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 1;
					buf1 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 2;
					buf2 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 3;
					buf3 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 4;
					buf4 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 5;
					buf5 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 6;
					buf6 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 7;
					buf7 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 8;
					buf8 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					i = i1 + 9;
					buf9 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]);
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {

						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		}


		fprintf(fp, "\n");

		// Heat Flux Y
		for (integer i1 = 0; i1 < maxelm; i1++) {
			if ((i1 + 10) < maxelm) {

				i = i1;
				buf0 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 1;
				buf1 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 2;
				buf2 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 3;
				buf3 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 4;
				buf4 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 5;
				buf5 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 6;
				buf6 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 7;
				buf7 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 8;
				buf8 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				i = i1 + 9;
				buf9 = signlog10(-t.prop[LAM][i] * Ty[0][i]);

				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

				i1 += 9;
			}
			else {
				i = i1;
				buf0 = signlog10(-t.prop[LAM][i] * Ty[0][i]);
				fprintf(fp, "%+.16f ", buf0);
				if (i1 % 10 == 0) fprintf(fp, "\n");
			}
		}

		if (bextendedprint) {
			for (integer i1 = 0; i1 < t.maxbound; i1++) {
				if ((i1 + 10) < t.maxbound) {

					i = i1;
					buf0 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);


					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop_b[LAM][i] * Ty[0][i + maxelm]);
					fprintf(fp, "%+.16f ", buf0);
					if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// Heat Flux Y
			for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
				if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {

					i = i1;
					buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 1;
					buf1 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 2;
					buf2 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 3;
					buf3 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 4;
					buf4 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 5;
					buf5 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 6;
					buf6 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 7;
					buf7 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 8;
					buf8 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					i = i1 + 9;
					buf9 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]);
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {

						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		}

		fprintf(fp, "\n");

		// Heat Flux Z
		for (integer i1 = 0; i1 < maxelm; i1++) {
			if ((i1 + 10) < maxelm) {

				i = i1;
				buf0 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 1;
				buf1 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 2;
				buf2 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 3;
				buf3 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 4;
				buf4 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 5;
				buf5 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 6;
				buf6 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 7;
				buf7 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 8;
				buf8 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				i = i1 + 9;
				buf9 = signlog10(-t.prop[LAM][i] * Tz[0][i]);

				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

				i1 += 9;
			}
			else {
				i = i1;
				buf0 = signlog10(-t.prop[LAM][i] * Tz[0][i]);
				fprintf(fp, "%+.16f ", buf0);
				if (i1 % 10 == 0) fprintf(fp, "\n");
			}
		}

		if (bextendedprint) {
			for (integer i1 = 0; i1 < t.maxbound; i1++) {
				if ((i1 + 10) < t.maxbound) {

					i = i1;
					buf0 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 1;
					buf1 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 2;
					buf2 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 3;
					buf3 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 4;
					buf4 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 5;
					buf5 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 6;
					buf6 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 7;
					buf7 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 8;
					buf8 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					i = i1 + 9;
					buf9 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);


					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-t.prop_b[LAM][i] * Tz[0][i + maxelm]);
					fprintf(fp, "%+.16f ", buf0);
					if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// Heat Flux Z
			for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
				if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {

					i = i1;
					buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 1;
					buf1 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 2;
					buf2 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 3;
					buf3 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 4;
					buf4 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 5;
					buf5 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 6;
					buf6 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 7;
					buf7 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 8;
					buf8 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					i = i1 + 9;
					buf9 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = signlog10(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]);
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {

						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 1;
						buf1 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 2;
						buf2 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 3;
						buf3 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 4;
						buf4 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 5;
						buf5 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 6;
						buf6 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 7;
						buf7 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 8;
						buf8 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						i = i1 + 9;
						buf9 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						fprintf(fp, "%+.16f ", buf0);
						if ((i1 + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		}

		fprintf(fp, "\n");


		// Mag Heat Flux
		for (integer i1 = 0; i1 < maxelm; i1++) {
			if ((i1 + 10) < maxelm) {
				i = i1;
				buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 1;
				buf1 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 2;
				buf2 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 3;
				buf3 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 4;
				buf4 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 5;
				buf5 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 6;
				buf6 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 7;
				buf7 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 8;
				buf8 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 9;
				buf9 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));


				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

				i1 += 9;
			}
			else {
				i = i1;
				buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				fprintf(fp, "%+.16f ", buf0);
				if (i1 % 10 == 0) fprintf(fp, "\n");
			}
		}

		if (bextendedprint) {
			for (integer i1 = 0; i1 < t.maxbound; i1++) {
				if ((i1 + 10) < t.maxbound) {
					i = i1;
					buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					integer i = i1;
					buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
					fprintf(fp, "%+.16f ", buf0);


					if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// Mag Heat Flux
			for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
				if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {
					i = i1;
					buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 1;
					buf1 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 2;
					buf2 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 3;
					buf3 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 4;
					buf4 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 5;
					buf5 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 6;
					buf6 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 7;
					buf7 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 8;
					buf8 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 9;
					buf9 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));


					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {
						i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						//fprintf(fp, "%+.16f ", -my_union[iunion_scan].t.prop_b[LAM][i] * Tx[i + maxelm]);
						fprintf(fp, "%+.16f ", buf0);


						if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		}


		fprintf(fp, "\n");


		// log10 Mag Heat Flux
		const doublereal eps_min = 2.0;
		const doublereal d_stub = 0.0;
		for (integer i1 = 0; i1 < maxelm; i1++) {
			if ((i1 + 10) < maxelm) {
				i = i1;
				buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 1;
				buf1 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 2;
				buf2 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 3;
				buf3 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 4;
				buf4 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 5;
				buf5 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 6;
				buf6 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 7;
				buf7 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 8;
				buf8 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				i = i1 + 9;
				buf9 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));

				if (buf0 > eps_min) {
					buf0 = log10(buf0);
				}
				else {
					buf0 = d_stub;
				}
				if (buf1 > eps_min) {
					buf1 = log10(buf1);
				}
				else {
					buf1 = d_stub;
				}
				if (buf2 > eps_min) {
					buf2 = log10(buf2);
				}
				else {
					buf2 = d_stub;
				}
				if (buf3 > eps_min) {
					buf3 = log10(buf3);
				}
				else {
					buf3 = d_stub;
				}
				if (buf4 > eps_min) {
					buf4 = log10(buf4);
				}
				else {
					buf4 = d_stub;
				}
				if (buf5 > eps_min) {
					buf5 = log10(buf5);
				}
				else {
					buf5 = d_stub;
				}
				if (buf6 > eps_min) {
					buf6 = log10(buf6);
				}
				else {
					buf6 = d_stub;
				}
				if (buf7 > eps_min) {
					buf7 = log10(buf7);
				}
				else {
					buf7 = d_stub;
				}
				if (buf8 > eps_min) {
					buf8 = log10(buf8);
				}
				else {
					buf8 = d_stub;
				}
				if (buf9 > eps_min) {
					buf9 = log10(buf9);
				}
				else {
					buf9 = d_stub;
				}

				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

				i1 += 9;
			}
			else {
				i = i1;
				buf0 = sqrt((-t.prop[LAM][i] * Tx[0][i])*(-t.prop[LAM][i] * Tx[0][i]) + (-t.prop[LAM][i] * Ty[0][i])*(-t.prop[LAM][i] * Ty[0][i]) + (-t.prop[LAM][i] * Tz[0][i])*(-t.prop[LAM][i] * Tz[0][i]));
				if (buf0 > eps_min) {
					buf0 = log10(buf0);
				}
				else {
					buf0 = d_stub;
				}
				fprintf(fp, "%+.16f ", buf0);
				if (i1 % 10 == 0) fprintf(fp, "\n");
			}
		}

		if (bextendedprint) {
			for (integer i1 = 0; i1 < t.maxbound; i1++) {
				if ((i1 + 10) < t.maxbound) {
					i = i1;
					buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 1;
					buf1 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 2;
					buf2 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 3;
					buf3 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 4;
					buf4 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 5;
					buf5 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 6;
					buf6 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 7;
					buf7 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 8;
					buf8 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));
					i = i1 + 9;
					buf9 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));

					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					if (buf1 > eps_min) {
						buf1 = log10(buf1);
					}
					else {
						buf1 = d_stub;
					}
					if (buf2 > eps_min) {
						buf2 = log10(buf2);
					}
					else {
						buf2 = d_stub;
					}
					if (buf3 > eps_min) {
						buf3 = log10(buf3);
					}
					else {
						buf3 = d_stub;
					}
					if (buf4 > eps_min) {
						buf4 = log10(buf4);
					}
					else {
						buf4 = d_stub;
					}
					if (buf5 > eps_min) {
						buf5 = log10(buf5);
					}
					else {
						buf5 = d_stub;
					}
					if (buf6 > eps_min) {
						buf6 = log10(buf6);
					}
					else {
						buf6 = d_stub;
					}
					if (buf7 > eps_min) {
						buf7 = log10(buf7);
					}
					else {
						buf7 = d_stub;
					}
					if (buf8 > eps_min) {
						buf8 = log10(buf8);
					}
					else {
						buf8 = d_stub;
					}
					if (buf9 > eps_min) {
						buf9 = log10(buf9);
					}
					else {
						buf9 = d_stub;
					}

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					integer i = i1;
					buf0 = sqrt((-t.prop_b[LAM][i] * Tx[0][i + maxelm])*(-t.prop_b[LAM][i] * Tx[0][i + maxelm]) + (-t.prop_b[LAM][i] * Ty[0][i + maxelm])*(-t.prop_b[LAM][i] * Ty[0][i + maxelm]) + (-t.prop_b[LAM][i] * Tz[0][i + maxelm])*(-t.prop_b[LAM][i] * Tz[0][i + maxelm]));

					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}

					//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[0][i +maxelm]);
					fprintf(fp, "%+.16f ", buf0);


					if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
		}

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// log10 Mag Heat Flux
			const doublereal eps_min = 2.0;
			const doublereal d_stub = 0.0;
			for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxelm; i1++) {
				if ((i1 + 10) < my_union[iunion_scan].t.maxelm) {
					i = i1;
					buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 1;
					buf1 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 2;
					buf2 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 3;
					buf3 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 4;
					buf4 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 5;
					buf5 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 6;
					buf6 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 7;
					buf7 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 8;
					buf8 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					i = i1 + 9;
					buf9 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));

					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					if (buf1 > eps_min) {
						buf1 = log10(buf1);
					}
					else {
						buf1 = d_stub;
					}
					if (buf2 > eps_min) {
						buf2 = log10(buf2);
					}
					else {
						buf2 = d_stub;
					}
					if (buf3 > eps_min) {
						buf3 = log10(buf3);
					}
					else {
						buf3 = d_stub;
					}
					if (buf4 > eps_min) {
						buf4 = log10(buf4);
					}
					else {
						buf4 = d_stub;
					}
					if (buf5 > eps_min) {
						buf5 = log10(buf5);
					}
					else {
						buf5 = d_stub;
					}
					if (buf6 > eps_min) {
						buf6 = log10(buf6);
					}
					else {
						buf6 = d_stub;
					}
					if (buf7 > eps_min) {
						buf7 = log10(buf7);
					}
					else {
						buf7 = d_stub;
					}
					if (buf8 > eps_min) {
						buf8 = log10(buf8);
					}
					else {
						buf8 = d_stub;
					}
					if (buf9 > eps_min) {
						buf9 = log10(buf9);
					}
					else {
						buf9 = d_stub;
					}

					fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

					i1 += 9;
				}
				else {
					i = i1;
					buf0 = sqrt((-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tx[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Ty[iunion_scan + 1][i]) + (-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i])*(-my_union[iunion_scan].t.prop[LAM][i] * Tz[iunion_scan + 1][i]));
					if (buf0 > eps_min) {
						buf0 = log10(buf0);
					}
					else {
						buf0 = d_stub;
					}
					fprintf(fp, "%+.16f ", buf0);
					if (i1 % 10 == 0) fprintf(fp, "\n");
				}
			}

			if (bextendedprint) {
				for (integer i1 = 0; i1 < my_union[iunion_scan].t.maxbound; i1++) {
					if ((i1 + 10) < my_union[iunion_scan].t.maxbound) {
						i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 1;
						buf1 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 2;
						buf2 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 3;
						buf3 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 4;
						buf4 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 5;
						buf5 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 6;
						buf6 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 7;
						buf7 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 8;
						buf8 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));
						i = i1 + 9;
						buf9 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						integer i = i1;
						buf0 = sqrt((-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Ty[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]) + (-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm])*(-my_union[iunion_scan].t.prop_b[LAM][i] * Tz[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}

						//fprintf(fp, "%+.16f ", -my_union[iunion_scan].t.prop_b[LAM][i] * Tx[iunion_scan + 1][i + my_union[iunion_scan].t.maxelm]);
						fprintf(fp, "%+.16f ", buf0);


						if ((i1 + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		}

		fprintf(fp, "\n");


		}


		// Освобождение оперативной памяти.

		for (integer iuscan = 0; iuscan <= lu; iuscan++) {
			if (Tx[iuscan] != nullptr) {
				delete[] Tx[iuscan];
			}
			if (Ty[iuscan] != nullptr) {
				delete[] Ty[iuscan];
			}
			if (Tz[iuscan] != nullptr) {
				delete[] Tz[iuscan];
			}
		}
		if (Tx != nullptr) {
			delete[] Tx;
		}
		if (Ty != nullptr) {
			delete[] Ty;
		}
		if (Tz != nullptr) {
			delete[] Tz;
		}

		fprintf(fp, "\n");

		for (integer j_6 = 0; j_6 < 4; j_6++) {

			// запись полной деформации
			if (lite_export) {
				for (i = 0; i < maxelm; i++) {


					if (ionly_solid_visible == 1) {
						fprintf(fp, "%e ", total_deformation_shadow[0][j_6][i]);
						//printf("%e \n", total_deformation_shadow[j_6][i]);
					}
					else {
						fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						//printf("%e \n", t.total_deformation[j_6][i]);
					}
					if (i % 10 == 0) {
						fprintf(fp, "\n");
						//getchar();
					}
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[0][j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
			else {
				for (i = 0; i < maxelm; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%e ", total_deformation_shadow[0][j_6][i]);
					}
					else {
						fprintf(fp, "%e ", t.total_deformation[j_6][i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[0][j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// запись полной деформации
				if (lite_export) {
					for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {


						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[iunion_scan + 1][j_6][i]);
							//printf("%e \n", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", my_union[iunion_scan].t.total_deformation[j_6][i]);
							//printf("%e \n", t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) {
							fprintf(fp, "\n");
							//getchar();
						}
					}

					if (bextendedprint) {
						for (i = my_union[iunion_scan].t.maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
							if (ionly_solid_visible == 1) {
								fprintf(fp, "%e ", total_deformation_shadow[iunion_scan + 1][j_6][i]);
							}
							else {
								fprintf(fp, "%e ", my_union[iunion_scan].t.total_deformation[j_6][i]);
							}
							if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[iunion_scan + 1][j_6][i]);
						}
						else {
							fprintf(fp, "%e ", my_union[iunion_scan].t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (bextendedprint) {
						for (i = maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
							if (ionly_solid_visible == 1) {
								fprintf(fp, "%e ", total_deformation_shadow[iunion_scan + 1][j_6][i]);
							}
							else {
								fprintf(fp, "%e ", my_union[iunion_scan].t.total_deformation[j_6][i]);
							}
							if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
			}
		

			fprintf(fp, "\n");
        }
		
		

		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					/*
					TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
					*/
					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					/*
					in_model_temp(p1, ib1, b, lb);
					in_model_temp(p2, ib2, b, lb);
					in_model_temp(p3, ib3, b, lb);
					in_model_temp(p4, ib4, b, lb);
					in_model_temp(p5, ib5, b, lb);
					in_model_temp(p6, ib6, b, lb);
					in_model_temp(p7, ib7, b, lb);
					in_model_temp(p8, ib8, b, lb);
					*/

					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];


					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
						integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
						integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
						integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


						integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
						integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
						integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
						integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



						integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
						integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
						integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
						integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;


						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
						pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
						pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);

						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
					}
					else {
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);

						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}

			// TODO
			// Организовать счётчик считающий количество ячеек.
			integer iadd = t.maxelm;
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				if (iunion_scan > 0) {
					iadd += my_union[iunion_scan-1].t.maxelm;
				}

				// запись информации о разностной сетке
				for (i = 0; i < my_union[iunion_scan].t.database.ncell; i++) {
					if (bsolid_static_only) {
						//printf("Only solid ok\n");
						//getchar();

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = my_union[iunion_scan].t.database.nvtxcell[0][i] - 1;
						inode2 = my_union[iunion_scan].t.database.nvtxcell[1][i] - 1;
						inode3 = my_union[iunion_scan].t.database.nvtxcell[2][i] - 1;
						inode4 = my_union[iunion_scan].t.database.nvtxcell[3][i] - 1;
						inode5 = my_union[iunion_scan].t.database.nvtxcell[4][i] - 1;
						inode6 = my_union[iunion_scan].t.database.nvtxcell[5][i] - 1;
						inode7 = my_union[iunion_scan].t.database.nvtxcell[6][i] - 1;
						inode8 = my_union[iunion_scan].t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						/*
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);
						*/

						ib1 = my_union[iunion_scan].t.whot_is_block[inode1];
						ib2 = my_union[iunion_scan].t.whot_is_block[inode2];
						ib3 = my_union[iunion_scan].t.whot_is_block[inode3];
						ib4 = my_union[iunion_scan].t.whot_is_block[inode4];
						ib5 = my_union[iunion_scan].t.whot_is_block[inode5];
						ib6 = my_union[iunion_scan].t.whot_is_block[inode6];
						ib7 = my_union[iunion_scan].t.whot_is_block[inode7];
						ib8 = my_union[iunion_scan].t.whot_is_block[inode8];


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// Визуализация твёрдого тела и только как и раньше.
							//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i], iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и только как и раньше.
							//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i], iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif


						}

					}
					else {
						//printf("fluid plot\n");
						//getchar();

						if ((ionly_solid_visible == 1) && (flow_interior > 0))
						{
							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = my_union[iunion_scan].t.database.nvtxcell[0][i] - 1;
							inode2 = my_union[iunion_scan].t.database.nvtxcell[1][i] - 1;
							inode3 = my_union[iunion_scan].t.database.nvtxcell[2][i] - 1;
							inode4 = my_union[iunion_scan].t.database.nvtxcell[3][i] - 1;
							inode5 = my_union[iunion_scan].t.database.nvtxcell[4][i] - 1;
							inode6 = my_union[iunion_scan].t.database.nvtxcell[5][i] - 1;
							inode7 = my_union[iunion_scan].t.database.nvtxcell[6][i] - 1;
							inode8 = my_union[iunion_scan].t.database.nvtxcell[7][i] - 1;

							integer inode2W = my_union[iunion_scan].t.sosedi[WSIDE][inode1].iNODE1;
							integer inode3W = my_union[iunion_scan].t.sosedi[WSIDE][inode4].iNODE1;
							integer inode6W = my_union[iunion_scan].t.sosedi[WSIDE][inode5].iNODE1;
							integer inode7W = my_union[iunion_scan].t.sosedi[WSIDE][inode8].iNODE1;


							integer inode5B = my_union[iunion_scan].t.sosedi[BSIDE][inode1].iNODE1;
							integer inode6B = my_union[iunion_scan].t.sosedi[BSIDE][inode2].iNODE1;
							integer inode7B = my_union[iunion_scan].t.sosedi[BSIDE][inode3].iNODE1;
							integer inode8B = my_union[iunion_scan].t.sosedi[BSIDE][inode4].iNODE1;



							integer inode3S = my_union[iunion_scan].t.sosedi[SSIDE][inode2].iNODE1;
							integer inode4S = my_union[iunion_scan].t.sosedi[SSIDE][inode1].iNODE1;
							integer inode7S = my_union[iunion_scan].t.sosedi[SSIDE][inode6].iNODE1;
							integer inode8S = my_union[iunion_scan].t.sosedi[SSIDE][inode5].iNODE1;


							TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
							center_cord3D(inode1, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p1, 100);
							center_cord3D(inode2, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p2, 100);
							center_cord3D(inode3, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p3, 100);
							center_cord3D(inode4, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p4, 100);
							center_cord3D(inode5, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p5, 100);
							center_cord3D(inode6, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p6, 100);
							center_cord3D(inode7, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p7, 100);
							center_cord3D(inode8, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p8, 100);
							pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
							pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
							pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
							in_model_temp(p1, ib1, b, lb);
							in_model_temp(p2, ib2, b, lb);
							in_model_temp(p3, ib3, b, lb);
							in_model_temp(p4, ib4, b, lb);
							in_model_temp(p5, ib5, b, lb);
							in_model_temp(p6, ib6, b, lb);
							in_model_temp(p7, ib7, b, lb);
							in_model_temp(p8, ib8, b, lb);

							// визуализация только твёрдого тела.
							// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
							if (((my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) &&
								(my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))
							{

								if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
									fprintf(fp, "%lld %lld %lld %lld ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%lld %lld %lld %lld\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif

								}
							}
							else if (((inode5B >= 0) && (inode5B < my_union[iunion_scan].t.maxelm) && (inode6B >= 0) && (inode6B < my_union[iunion_scan].t.maxelm) && (inode7B >= 0) && (inode7B < my_union[iunion_scan].t.maxelm) && (inode8B >= 0) && (inode8B < my_union[iunion_scan].t.maxelm) && (my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) && (!((my_union[iunion_scan].t.ptr[1][inode5B] == -1) && (my_union[iunion_scan].t.ptr[1][inode6B] == -1) && (my_union[iunion_scan].t.ptr[1][inode7B] == -1) && (my_union[iunion_scan].t.ptr[1][inode8B] == -1))) && (!((my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))))
							{
								if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
									fprintf(fp, "%lld %lld %lld %lld ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%lld %lld %lld %lld\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif

								}
							}
							else if (((inode2W >= 0) && (inode2W < my_union[iunion_scan].t.maxelm) && (inode3W >= 0) && (inode3W < my_union[iunion_scan].t.maxelm) && (inode6W >= 0) && (inode6W < my_union[iunion_scan].t.maxelm) && (inode7W >= 0) && (inode7W < my_union[iunion_scan].t.maxelm) && (my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) && (my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1) && (!((my_union[iunion_scan].t.ptr[1][inode2W] == -1) && (my_union[iunion_scan].t.ptr[1][inode3W] == -1) && (my_union[iunion_scan].t.ptr[1][inode6W] == -1) && (my_union[iunion_scan].t.ptr[1][inode7W] == -1))) && (!((my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1)))))
							{
								if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
									fprintf(fp, "%lld %lld %lld %lld ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%lld %lld %lld %lld\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif

								}
							}
							else if (((inode3S >= 0) && (inode3S < my_union[iunion_scan].t.maxelm) && (inode4S >= 0) && (inode4S < my_union[iunion_scan].t.maxelm) && (inode7S >= 0) && (inode7S < my_union[iunion_scan].t.maxelm) && (inode8S >= 0) && (inode8S < my_union[iunion_scan].t.maxelm) && (my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (!((my_union[iunion_scan].t.ptr[1][inode3S] == -1) && (my_union[iunion_scan].t.ptr[1][inode4S] == -1) && (my_union[iunion_scan].t.ptr[1][inode7S] == -1) && (my_union[iunion_scan].t.ptr[1][inode8S] == -1))) && (!((my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))))
							{
								if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
									fprintf(fp, "%lld %lld %lld %lld ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%lld %lld %lld %lld\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif

								}
							}
						}
						else {
							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = my_union[iunion_scan].t.database.nvtxcell[0][i] - 1;
							inode2 = my_union[iunion_scan].t.database.nvtxcell[1][i] - 1;
							inode3 = my_union[iunion_scan].t.database.nvtxcell[2][i] - 1;
							inode4 = my_union[iunion_scan].t.database.nvtxcell[3][i] - 1;
							inode5 = my_union[iunion_scan].t.database.nvtxcell[4][i] - 1;
							inode6 = my_union[iunion_scan].t.database.nvtxcell[5][i] - 1;
							inode7 = my_union[iunion_scan].t.database.nvtxcell[6][i] - 1;
							inode8 = my_union[iunion_scan].t.database.nvtxcell[7][i] - 1;

							TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
							//TOCHKA pall;
							center_cord3D(inode1, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p1, 100);
							center_cord3D(inode2, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p2, 100);
							center_cord3D(inode3, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p3, 100);
							center_cord3D(inode4, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p4, 100);
							center_cord3D(inode5, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p5, 100);
							center_cord3D(inode6, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p6, 100);
							center_cord3D(inode7, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p7, 100);
							center_cord3D(inode8, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa, p8, 100);

							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
							in_model_temp(p1, ib1, b, lb);
							in_model_temp(p2, ib2, b, lb);
							in_model_temp(p3, ib3, b, lb);
							in_model_temp(p4, ib4, b, lb);
							in_model_temp(p5, ib5, b, lb);
							in_model_temp(p6, ib6, b, lb);
							in_model_temp(p7, ib7, b, lb);
							in_model_temp(p8, ib8, b, lb);

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								// Визуализация твёрдого тела и жидкости.
								// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd+my_union[iunion_scan].t.database.nvtxcell[1][i], iadd+my_union[iunion_scan].t.database.nvtxcell[2][i], iadd+my_union[iunion_scan].t.database.nvtxcell[3][i], iadd+my_union[iunion_scan].t.database.nvtxcell[4][i], iadd+my_union[iunion_scan].t.database.nvtxcell[5][i], iadd+my_union[iunion_scan].t.database.nvtxcell[6][i], iadd+my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
								// Визуализация твёрдого тела и жидкости.
								// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								fprintf(fp, "%d %d %d %d %d %d %d %d\n", iadd+my_union[iunion_scan].t.database.nvtxcell[0][i], iadd+my_union[iunion_scan].t.database.nvtxcell[1][i], iadd+my_union[iunion_scan].t.database.nvtxcell[2][i], iadd+my_union[iunion_scan].t.database.nvtxcell[3][i], iadd+my_union[iunion_scan].t.database.nvtxcell[4][i], iadd+my_union[iunion_scan].t.database.nvtxcell[5][i], iadd+my_union[iunion_scan].t.database.nvtxcell[6][i], iadd+my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif


							}

						}
					}
				}
			}

		}
		else {
#ifdef MINGW_COMPILLER
		err = 0;
		fp1=fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
		if (fp1 == nullptr) err = 1;
#else
		err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	for (integer iunion_scan = 0; iunion_scan <= lu; iunion_scan++) {
		if (temp_shadow[iunion_scan] != nullptr) {
			delete[] temp_shadow[iunion_scan];
			temp_shadow[iunion_scan] = nullptr;
		}
	}

	if (temp_shadow != nullptr) {
		delete[] temp_shadow;
		temp_shadow = nullptr;
	}
	//total_deformation

	if (ionly_solid_visible == 1) {
		if (total_deformation_shadow != nullptr) {

			for (integer iunion_scan = 0; iunion_scan <= lu; iunion_scan++) {
				for (integer j_6 = 0; j_6 < 4; j_6++) {
					if (total_deformation_shadow[iunion_scan][j_6] != nullptr) {
						delete[] total_deformation_shadow[iunion_scan][j_6];
						total_deformation_shadow[iunion_scan][j_6] = nullptr;
					}
				}
			}

			for (integer iunion_scan = 0; iunion_scan <= lu; iunion_scan++) {
				delete[] total_deformation_shadow[iunion_scan];
				total_deformation_shadow[iunion_scan] = nullptr;
			}

			delete[] total_deformation_shadow;
			total_deformation_shadow = nullptr;
		}
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);



} // exporttecplotxy360T_3D_part2_assembles

// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2. Для анимации. Если inumbercadr==0 то первый кадр.
void exporttecplotxy360T_3D_part2_ianimation_series( integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, integer ikey,
	integer inumbercadr, doublereal time, BLOCK* b)
{
	const bool lite_export = true;
	// 16 знаков после запятой сохранять никому ненужно,
	// вполне достаточно шести знаков.

#if doubleintprecision == 1
	printf("ionly_solid_visible =%lld\n", ionly_solid_visible);
#else
	printf("ionly_solid_visible =%d\n", ionly_solid_visible);
#endif


	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal* temp_shadow = nullptr;
	if (ionly_solid_visible == 1) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == 1) {
		total_deformation_shadow = new doublereal*[4];
		for (integer j_1 = 0; j_1 < 4; j_1++) {
			total_deformation_shadow[j_1] = new doublereal[t.maxelm + t.maxbound];
			for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
				total_deformation_shadow[j_1][i_1] = t.total_deformation[j_1][i_1];
			}
		}
	}


	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	err = 0;
	if (inumbercadr == 0) {
		switch (ikey) {
		case 0: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		case 1: fp = fopen64("ALICEFLOW0_27_temp_animation.PLT", "wb");  break; // То что нужно для отчёта Алексею.
		default: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		}
}
	else {
		switch (ikey) {
		case 0: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		case 1: fp = fopen64("ALICEFLOW0_27_temp_animation.PLT", "ab");  break; // То что нужно для отчёта Алексею.
		default: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		}
	}
	if (fp==NULL) err = 1;
#else
	if (inumbercadr == 0) {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_animation.PLT", "wb");  break; // То что нужно для отчёта Алексею.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		}
	}
	else {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_animation.PLT", "ab");  break; // То что нужно для отчёта Алексею.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		}
	}
#endif
	
	
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2_ianimation_series in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1 = fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1==NULL) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		if (bOk)
		{


			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count > 0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i = 0; i < flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // считаем что температура всегда активна
				}
			}

			if (ivarexport == 1) {

				integer ncell_shadow = 0;
				for (i = 0; i < ncell; i++) {

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					// Использование быстродействующей хеш таблицы whot_is_block значительно
					// быстрее и не приводит к каким бы то ни было отличиям от прямого метда.
					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];

					//in_model_temp(p1, ib1, b, lb);
					//in_model_temp(p2, ib2, b, lb);
					//in_model_temp(p3, ib3, b, lb);
					//in_model_temp(p4, ib4, b, lb);
					//in_model_temp(p5, ib5, b, lb);
					//in_model_temp(p6, ib6, b, lb);
					//in_model_temp(p7, ib7, b, lb);
					//in_model_temp(p8, ib8, b, lb);

					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
						ncell_shadow++;
					}
				}
				if (inumbercadr == 0) {
					// запись заголовка
					fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

					// запись имён переменных
					//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
					fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // печатается только поле температур
					

				}
#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm, ncell_shadow);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n",time,   maxelm, ncell_shadow);
				}
#endif

				// Уникальная подпись кадра.
				//fprintf(fp, "$!ATTACHTEXT\n");
				//fprintf(fp, "POSITIONCOORDSYS = GRID\n ANCHORPOS\n {\n X=%e\n Y=%e\ Z=%e\n }\n",b[0].g.xS+0.05*(b[0].g.xE-b[0].g.xS), b[0].g.yS + 0.05*(b[0].g.yE - b[0].g.yS), b[0].g.zS + 0.05*(b[0].g.zE - b[0].g.zS));
				//fprintf(fp, "TEXTSHAPE\n {\n ISBOLD = NO\n SIZEUNITS = GRID\n HEIGHT = %e\n }\n",0.07*sqrt((b[0].g.xE - b[0].g.xS)*(b[0].g.xE - b[0].g.xS)+(b[0].g.yE - b[0].g.yS)*(b[0].g.yE - b[0].g.yS)+(b[0].g.zE - b[0].g.zS)*(b[0].g.zE - b[0].g.zS)));
				//fprintf(fp, "BOX\n{\nBOXTYPE = FILLED\n	MARGIN = 30\n COLOR = WHITE\n }\n");
				//fprintf(fp, "TEXT='time=%e s'\n\n", time);



				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			else if (ivarexport == 3) {
				if (inumbercadr == 0) {
					// запись заголовка
					fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");
				}

				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ничего не делаем
					ncell_shadow = 0;
					for (i = 0; i < ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						/*
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);
						*/

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];


						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
							ncell_shadow++;
						}
					}
					//ncell_shadow = ncell;
				}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						ncell_shadow = 0;
						for (i = 0; i < t.database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = t.database.nvtxcell[0][i] - 1;
							inode2 = t.database.nvtxcell[1][i] - 1;
							inode3 = t.database.nvtxcell[2][i] - 1;
							inode4 = t.database.nvtxcell[3][i] - 1;
							inode5 = t.database.nvtxcell[4][i] - 1;
							inode6 = t.database.nvtxcell[5][i] - 1;
							inode7 = t.database.nvtxcell[6][i] - 1;
							inode8 = t.database.nvtxcell[7][i] - 1;

							integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
							integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
							integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
							integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


							integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
							integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
							integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
							integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;


							integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
							integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
							integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
							integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;

							//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
							//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
							//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
							//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
							//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
							//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
							//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
							//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
							//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
							//pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
							//pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
							//pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

							integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;

							//in_model_temp(p1, ib1, b, lb);
							//in_model_temp(p2, ib2, b, lb);
							//in_model_temp(p3, ib3, b, lb);
							//in_model_temp(p4, ib4, b, lb);
							//in_model_temp(p5, ib5, b, lb);
							//in_model_temp(p6, ib6, b, lb);
							//in_model_temp(p7, ib7, b, lb);
							//in_model_temp(p8, ib8, b, lb);

							ib1 = t.whot_is_block[inode1];
							ib2 = t.whot_is_block[inode2];
							ib3 = t.whot_is_block[inode3];
							ib4 = t.whot_is_block[inode4];
							ib5 = t.whot_is_block[inode5];
							ib6 = t.whot_is_block[inode6];
							ib7 = t.whot_is_block[inode7];
							ib8 = t.whot_is_block[inode8];

							// визуализация только твёрдого тела.
							// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
							if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
								&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
							{
								if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
									ncell_shadow++;
								}

							}
							else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
							{

								if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {
									ncell_shadow++;
									temp_shadow[inode5] = t.potent[inode1];
									temp_shadow[inode6] = t.potent[inode2];
									temp_shadow[inode7] = t.potent[inode3];
									temp_shadow[inode8] = t.potent[inode4];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[j_4][inode5] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode2];
										total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode3];
										total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode4];
									}
								}
							}
							else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
							{

								if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

									ncell_shadow++;
									temp_shadow[inode2] = t.potent[inode1];
									temp_shadow[inode3] = t.potent[inode4];
									temp_shadow[inode6] = t.potent[inode5];
									temp_shadow[inode7] = t.potent[inode8];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[j_4][inode2] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode4];
										total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode5];
										total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode8];
									}
								}
							}
							else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
							{

								if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

									ncell_shadow++;
									temp_shadow[inode4] = t.potent[inode1];
									temp_shadow[inode3] = t.potent[inode2];
									temp_shadow[inode8] = t.potent[inode5];
									temp_shadow[inode7] = t.potent[inode6];
									// total_deformation
									for (integer j_4 = 0; j_4 < 4; j_4++) {
										total_deformation_shadow[j_4][inode4] = t.total_deformation[j_4][inode1];
										total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode2];
										total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode5];
										total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode6];
									}
								}
							}
						}
						if (ncell_shadow == 0) {
							// У нас совсем нету твёрдого тела, поэтому мы показываем только  жидкость.
							ncell_shadow = ncell;
							ionly_solid_visible = 0;
						}
					}
				}
				if (inumbercadr == 0) {
					// Полный набор искомых величин и теплопередача и гидродинамика:
					if (bextendedprint) {
						//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
						fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
						
					}
					else {
						//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
						fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
						
					}
				}
#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm, ncell_shadow);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm, ncell_shadow);
				}
#endif



				
					if (bvery_big_memory) {
						if (lite_export) {
							// запись x
							for (i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись y
							for (i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись z
							for (i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						else {
							// запись x
							for (i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись y
							for (i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись z
							for (i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
					else {
						while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					}
				}
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		

		// запись второй части

		// Запись поля температур производится всегда.


		// запись температуры
		
			if (lite_export) {
				for (i = 0; i < maxelm; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.3f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.3f ", t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%+.3f ", temp_shadow[i]);
						}
						else {
							fprintf(fp, "%+.3f ", t.potent[i]);
						}
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
			else {
				for (i = 0; i < maxelm; i++) {
					if (ionly_solid_visible == 1) {
						fprintf(fp, "%+.16f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.16f ", t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%+.16f ", temp_shadow[i]);
						}
						else {
							fprintf(fp, "%+.16f ", t.potent[i]);
						}
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
		
			

		    fprintf(fp, "\n");
	    

		
			// Lam
			if (lite_export) {
				for (i = 0; i < maxelm; i++) {
					fprintf(fp, "%+.4f ", t.prop[LAM][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = 0; i < t.maxbound; i++) {
						fprintf(fp, "%+.4f ", t.prop_b[LAM][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
			else {
				for (i = 0; i < maxelm; i++) {
					fprintf(fp, "%+.16f ", t.prop[LAM][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = 0; i < t.maxbound; i++) {
						fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}


			fprintf(fp, "\n");
		

		
			// Запись гидродинамических величин если необходимо:
			if (ivarexport == 3) {
				// Speed
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						doublereal svx = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
						doublereal svy = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
						doublereal svz = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}


				if (bextendedprint) {
					// Speed
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						doublereal svx = f[idfluid].potent[VX][i + maxelm] * f[idfluid].potent[VX][i + maxelm];
						doublereal svy = f[idfluid].potent[VY][i + maxelm] * f[idfluid].potent[VY][i + maxelm];
						doublereal svz = f[idfluid].potent[VZ][i + maxelm] * f[idfluid].potent[VZ][i + maxelm];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}


				fprintf(fp, "\n");

				// Pressure
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Pressure
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// PAM
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// PAM
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// VX
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VX
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i + maxelm]); // VX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// VY
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VY
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i + maxelm]); // VY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// VZ
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VZ
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i + maxelm]); // VZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");


				// Rho
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
					}
					else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Rho
					for (i = 0; i < f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// Mu
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Mu
					for (i = 0; i < f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[0].prop_b[MU][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// Mut // Турбулентная динамическая вязкость
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);

						//fprintf(fp, "%+.16f ", 1.0*f[t.ptr[1][i]].icolor_different_fluid_domain[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// MUT
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						//fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
						fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm] / f[0].prop_b[MU][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// для отладки график распределения нумерации контрольных объёмов.
				// или распределение расстояния до стенки Distance_Wall.
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", doublereal(i));
						if ((f[t.ptr[1][i]].iflowregime == ZEROEQMOD) ||
							(f[t.ptr[1][i]].iflowregime == SMAGORINSKY) ||
							(f[t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES) ||
							(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) ||
							(f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
							fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
						}
						else fprintf(fp, "%+.16f ", 0.0);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Distance_Wall.
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						if ((f[0].iflowregime == ZEROEQMOD) ||
							(f[0].iflowregime == SMAGORINSKY)||
							(f[0].iflowregime == RANS_SPALART_ALLMARES) ||
							(f[0].iflowregime == RANS_MENTER_SST) ||
							(f[0].iflowregime == RANS_STANDART_K_EPS)) {
							fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
						}
						else fprintf(fp, "%+.16f ", 0.0);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");



				// Curl // Завихрённость - модуль ротора скорости.
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Curl
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// Частные производные от компонент скорости !!!.

				// GRADXVX
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVX
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADYVX
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVX
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADZVX
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVX
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADXVY
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVY
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADYVY
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVY
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADZVY
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVY
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADXVZ
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVZ
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADYVZ
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVZ
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADZVZ
				for (i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVZ
					integer idfluid = 0;
					for (i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i + maxelm]); // GRADZVZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

			}
	

		doublereal *Tx = nullptr;
		doublereal *Ty = nullptr;
		doublereal *Tz = nullptr;
		Tx = new doublereal[t.maxelm + t.maxbound];
		Ty = new doublereal[t.maxelm + t.maxbound];
		Tz = new doublereal[t.maxelm + t.maxbound];

		// инициализация нулём.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// нахождение градиентов.
		for (i = 0; i<t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, false,
				t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
		}

		for (i = 0; i<t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, true,
				t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
		}


		// Сохранение в файл.
		
			if (lite_export) {
				doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

				// Heat Flux X
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop[LAM][i] * Tx[i]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop[LAM][i] * Tx[i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tx[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Y
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop[LAM][i] * Ty[i]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop[LAM][i] * Ty[i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Ty[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Z
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 1;
						buf1 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 2;
						buf2 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 3;
						buf3 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 4;
						buf4 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 5;
						buf5 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 6;
						buf6 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 7;
						buf7 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 8;
						buf8 = signlog10(-t.prop[LAM][i] * Tz[i]);
						i = i1 + 9;
						buf9 = signlog10(-t.prop[LAM][i] * Tz[i]);

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 1;
							buf1 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 2;
							buf2 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 3;
							buf3 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 4;
							buf4 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 5;
							buf5 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 6;
							buf6 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 7;
							buf7 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 8;
							buf8 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							i = i1 + 9;
							buf9 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);


							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = signlog10(-t.prop_b[LAM][i] * Tz[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// Mag Heat Flux
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));


						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// log10 Mag Heat Flux
				const doublereal eps_min = 2.0;
				const doublereal d_stub = 0.0;
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}
							if (buf1 > eps_min) {
								buf1 = log10(buf1);
							}
							else {
								buf1 = d_stub;
							}
							if (buf2 > eps_min) {
								buf2 = log10(buf2);
							}
							else {
								buf2 = d_stub;
							}
							if (buf3 > eps_min) {
								buf3 = log10(buf3);
							}
							else {
								buf3 = d_stub;
							}
							if (buf4 > eps_min) {
								buf4 = log10(buf4);
							}
							else {
								buf4 = d_stub;
							}
							if (buf5 > eps_min) {
								buf5 = log10(buf5);
							}
							else {
								buf5 = d_stub;
							}
							if (buf6 > eps_min) {
								buf6 = log10(buf6);
							}
							else {
								buf6 = d_stub;
							}
							if (buf7 > eps_min) {
								buf7 = log10(buf7);
							}
							else {
								buf7 = d_stub;
							}
							if (buf8 > eps_min) {
								buf8 = log10(buf8);
							}
							else {
								buf8 = d_stub;
							}
							if (buf9 > eps_min) {
								buf9 = log10(buf9);
							}
							else {
								buf9 = d_stub;
							}

							fprintf(fp, "%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}

							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.6f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

			}
			else {
				doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

				// Heat Flux X
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 1;
						buf1 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 2;
						buf2 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 3;
						buf3 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 4;
						buf4 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 5;
						buf5 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 6;
						buf6 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 7;
						buf7 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 8;
						buf8 = -t.prop[LAM][i] * Tx[i];
						i = i1 + 9;
						buf9 = -t.prop[LAM][i] * Tx[i];

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop[LAM][i] * Tx[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 1;
							buf1 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 2;
							buf2 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 3;
							buf3 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 4;
							buf4 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 5;
							buf5 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 6;
							buf6 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 7;
							buf7 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 8;
							buf8 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							i = i1 + 9;
							buf9 = -t.prop_b[LAM][i] * Tx[i + maxelm];


							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tx[i + maxelm];
							fprintf(fp, "%+.16f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Y
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 1;
						buf1 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 2;
						buf2 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 3;
						buf3 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 4;
						buf4 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 5;
						buf5 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 6;
						buf6 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 7;
						buf7 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 8;
						buf8 = -t.prop[LAM][i] * Ty[i];
						i = i1 + 9;
						buf9 = -t.prop[LAM][i] * Ty[i];

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop[LAM][i] * Ty[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 1;
							buf1 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 2;
							buf2 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 3;
							buf3 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 4;
							buf4 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 5;
							buf5 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 6;
							buf6 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 7;
							buf7 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 8;
							buf8 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							i = i1 + 9;
							buf9 = -t.prop_b[LAM][i] * Ty[i + maxelm];


							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = -t.prop_b[LAM][i] * Ty[i + maxelm];
							fprintf(fp, "%+.16f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");

				// Heat Flux Z
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						i = i1;
						buf0 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 1;
						buf1 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 2;
						buf2 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 3;
						buf3 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 4;
						buf4 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 5;
						buf5 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 6;
						buf6 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 7;
						buf7 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 8;
						buf8 = -t.prop[LAM][i] * Tz[i];
						i = i1 + 9;
						buf9 = -t.prop[LAM][i] * Tz[i];

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = -t.prop[LAM][i] * Tz[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 1;
							buf1 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 2;
							buf2 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 3;
							buf3 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 4;
							buf4 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 5;
							buf5 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 6;
							buf6 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 7;
							buf7 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 8;
							buf8 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							i = i1 + 9;
							buf9 = -t.prop_b[LAM][i] * Tz[i + maxelm];


							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							i = i1;
							buf0 = -t.prop_b[LAM][i] * Tz[i + maxelm];
							fprintf(fp, "%+.16f ", buf0);
							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// Mag Heat Flux
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));


						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.16f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


				// log10 Mag Heat Flux
				const doublereal eps_min = 2.0;
				const doublereal d_stub = 0.0;
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 1;
						buf1 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 2;
						buf2 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 3;
						buf3 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 4;
						buf4 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 5;
						buf5 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 6;
						buf6 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 7;
						buf7 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 8;
						buf8 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						i = i1 + 9;
						buf9 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));

						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						if (buf1 > eps_min) {
							buf1 = log10(buf1);
						}
						else {
							buf1 = d_stub;
						}
						if (buf2 > eps_min) {
							buf2 = log10(buf2);
						}
						else {
							buf2 = d_stub;
						}
						if (buf3 > eps_min) {
							buf3 = log10(buf3);
						}
						else {
							buf3 = d_stub;
						}
						if (buf4 > eps_min) {
							buf4 = log10(buf4);
						}
						else {
							buf4 = d_stub;
						}
						if (buf5 > eps_min) {
							buf5 = log10(buf5);
						}
						else {
							buf5 = d_stub;
						}
						if (buf6 > eps_min) {
							buf6 = log10(buf6);
						}
						else {
							buf6 = d_stub;
						}
						if (buf7 > eps_min) {
							buf7 = log10(buf7);
						}
						else {
							buf7 = d_stub;
						}
						if (buf8 > eps_min) {
							buf8 = log10(buf8);
						}
						else {
							buf8 = d_stub;
						}
						if (buf9 > eps_min) {
							buf9 = log10(buf9);
						}
						else {
							buf9 = d_stub;
						}

						fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

						i1 += 9;
					}
					else {
						i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						if (buf0 > eps_min) {
							buf0 = log10(buf0);
						}
						else {
							buf0 = d_stub;
						}
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 1;
							buf1 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 2;
							buf2 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 3;
							buf3 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 4;
							buf4 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 5;
							buf5 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 6;
							buf6 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 7;
							buf7 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 8;
							buf8 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));
							i = i1 + 9;
							buf9 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}
							if (buf1 > eps_min) {
								buf1 = log10(buf1);
							}
							else {
								buf1 = d_stub;
							}
							if (buf2 > eps_min) {
								buf2 = log10(buf2);
							}
							else {
								buf2 = d_stub;
							}
							if (buf3 > eps_min) {
								buf3 = log10(buf3);
							}
							else {
								buf3 = d_stub;
							}
							if (buf4 > eps_min) {
								buf4 = log10(buf4);
							}
							else {
								buf4 = d_stub;
							}
							if (buf5 > eps_min) {
								buf5 = log10(buf5);
							}
							else {
								buf5 = d_stub;
							}
							if (buf6 > eps_min) {
								buf6 = log10(buf6);
							}
							else {
								buf6 = d_stub;
							}
							if (buf7 > eps_min) {
								buf7 = log10(buf7);
							}
							else {
								buf7 = d_stub;
							}
							if (buf8 > eps_min) {
								buf8 = log10(buf8);
							}
							else {
								buf8 = d_stub;
							}
							if (buf9 > eps_min) {
								buf9 = log10(buf9);
							}
							else {
								buf9 = d_stub;
							}

							fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f \n", buf0, buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9);

							i1 += 9;
						}
						else {
							integer i = i1;
							buf0 = sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm]));

							if (buf0 > eps_min) {
								buf0 = log10(buf0);
							}
							else {
								buf0 = d_stub;
							}

							//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
							fprintf(fp, "%+.16f ", buf0);


							if ((i1 + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}

				fprintf(fp, "\n");


			
		}


		// Освобождение оперативной памяти.
		if (Tx != nullptr) {
			delete[] Tx;
		}
		if (Ty != nullptr) {
			delete[] Ty;
		}
		if (Tz != nullptr) {
			delete[] Tz;
		}

		
			fprintf(fp, "\n");

			for (integer j_6 = 0; j_6 < 4; j_6++) {

				// запись полной деформации
				if (lite_export) {
					for (i = 0; i < maxelm; i++) {


						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							//printf("%e \n", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							//printf("%e \n", t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) {
							fprintf(fp, "\n");
							//getchar();
						}
					}

					if (bextendedprint) {
						for (i = maxelm; i < maxelm + t.maxbound; i++) {
							if (ionly_solid_visible == 1) {
								fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							}
							else {
								fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							}
							if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					for (i = 0; i < maxelm; i++) {
						if (ionly_solid_visible == 1) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (bextendedprint) {
						for (i = maxelm; i < maxelm + t.maxbound; i++) {
							if (ionly_solid_visible == 1) {
								fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							}
							else {
								fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							}
							if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
						}
					}
				}


				fprintf(fp, "\n");

			}
		

		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					/*
					TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
					*/
					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					/*
					in_model_temp(p1, ib1, b, lb);
					in_model_temp(p2, ib2, b, lb);
					in_model_temp(p3, ib3, b, lb);
					in_model_temp(p4, ib4, b, lb);
					in_model_temp(p5, ib5, b, lb);
					in_model_temp(p6, ib6, b, lb);
					in_model_temp(p7, ib7, b, lb);
					in_model_temp(p8, ib8, b, lb);
					*/

					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];


					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
						integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
						integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
						integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


						integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
						integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
						integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
						integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



						integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
						integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
						integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
						integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;


						//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
						//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						//pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
						//pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
						//pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];

						//in_model_temp(p1, ib1, b, lb);
						//in_model_temp(p2, ib2, b, lb);
						//in_model_temp(p3, ib3, b, lb);
						//in_model_temp(p4, ib4, b, lb);
						//in_model_temp(p5, ib5, b, lb);
						//in_model_temp(p6, ib6, b, lb);
						//in_model_temp(p7, ib7, b, lb);
						//in_model_temp(p8, ib8, b, lb);

						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
					}
					else {
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						//center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						//center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						//center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						//center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						//center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						//center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						//center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						//center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;

						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];

						//in_model_temp(p1, ib1, b, lb);
						//in_model_temp(p2, ib2, b, lb);
						//in_model_temp(p3, ib3, b, lb);
						//in_model_temp(p4, ib4, b, lb);
						//in_model_temp(p5, ib5, b, lb);
						//in_model_temp(p6, ib6, b, lb);
						//in_model_temp(p7, ib7, b, lb);
						//in_model_temp(p8, ib8, b, lb);

						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}
		}
		else {
#ifdef MINGW_COMPILLER
		err = 0;
		fp1 = fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
		if (fp1 == nullptr) err = 1;
#else
		err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	if (temp_shadow != nullptr) {
		delete[] temp_shadow;
		temp_shadow = nullptr;
	}
	//total_deformation

	if (total_deformation_shadow != nullptr) {
		for (integer j_6 = 0; j_6 < 4; j_6++) {
			if (total_deformation_shadow[j_6] != nullptr) {
				delete[] total_deformation_shadow[j_6];
				total_deformation_shadow[j_6] = nullptr;
			}
		}
		delete[] total_deformation_shadow;
		total_deformation_shadow = nullptr;
	}


	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}

// специально для контроля над алгебраическим многосеточным методом.
// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2amg(TEMPER &t, doublereal* u, bool bextendedprint, integer imove)
{
	integer ianimate = 0;
	integer flow_interior_count = 1;
	// imove 0 или 1 для нумерации с нуля или с единицы.

	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // часть 1 или 3
	errno_t err=0;
#ifdef MINGW_COMPILLER
	err = 0;
	fp=fopen64("ALICEFLOW0_07_temp.PLT", "w");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
#endif
	
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.


	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2amg in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1 = fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1 == nullptr) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		if (bOk)
		{


			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // считаем что температура всегда активна
				}
			}

			if (ivarexport == 1) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// запись имён переменных
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
				fprintf(fp, "VARIABLES = x, y, z, u\n");  // печатается только передаваемое поле

				// запись информации о зонах
				if (bextendedprint) {
					printf("extended printeger ne predusmotreno\n");
					//getchar();
					system("PAUSE");
#if doubleintprecision == 1
					//	fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm + t.maxbound, ncell);
#else
					//	fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm + t.maxbound, ncell);
#endif

				}
				else {

#if doubleintprecision == 1
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm, t.database.ncell);
#else
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm, t.database.ncell);
#endif

						}

				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			else if (ivarexport == 3) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				/*
				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ничего не делаем
					ncell_shadow = ncell;
				}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						ncell_shadow = 0;
						for (i = 0; i < t.database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = t.database.nvtxcell[0][i] - 1;
							inode2 = t.database.nvtxcell[1][i] - 1;
							inode3 = t.database.nvtxcell[2][i] - 1;
							inode4 = t.database.nvtxcell[3][i] - 1;
							inode5 = t.database.nvtxcell[4][i] - 1;
							inode6 = t.database.nvtxcell[5][i] - 1;
							inode7 = t.database.nvtxcell[6][i] - 1;
							inode8 = t.database.nvtxcell[7][i] - 1;
							// визуализация только твёрдого тела.
							// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
							if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
								&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {
								ncell_shadow++;
							}
						}
						if (ncell_shadow == 0) {
							// У нас совсем нету твёрдого тела, поэтому мы показываем только  жидкость.
							ncell_shadow = ncell;
							ionly_solid_visible = 0;
						}
					}
				}
				*/
				// Полный набор искомых величин и теплопередача и гидродинамика:
				if (bextendedprint) {
					fprintf(fp, "\nVARIABLES = x, y, z, u\n");
				}
				else {
					fprintf(fp, "\nVARIABLES = x, y, z, u\n");
				}

				// запись информации о зонах
				if (bextendedprint) {
					printf("extended printeger ne predusmotreno\n");
					//getchar();
					system("PAUSE");

#if doubleintprecision == 1
					//fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm + t.maxbound, ncell_shadow);
#else
					//fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm + t.maxbound, ncell_shadow);
#endif

						}
				else {

#if doubleintprecision == 1
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm, t.database.ncell);
#else
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", t.database.maxelm, t.database.ncell);
#endif

					
				}
				if (bvery_big_memory) {
					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части

		// Запись поля температур производится всегда.

		/*
		// запись температуры
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.potent[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = maxelm; i<maxelm + t.maxbound; i++) {
				fprintf(fp, "%+.16f ", t.potent[i]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}


		fprintf(fp, "\n");

		// Lam
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.prop[LAM][i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");
		*/

		// Запись гидродинамических величин если необходимо:
		if (ivarexport == 3) {
			/*
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VX][i + maxelm] * f[idfluid].potent[VX][i + maxelm];
					doublereal svy = f[idfluid].potent[VY][i + maxelm] * f[idfluid].potent[VY][i + maxelm];
					doublereal svz = f[idfluid].potent[VZ][i + maxelm] * f[idfluid].potent[VZ][i + maxelm];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
			

			fprintf(fp, "\n");
			*/
			// Pressure
			doublereal avg_stat = 0.0;
			for (i = 0; i < t.database.maxelm; i++) {
				if (t.database.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
					fprintf(fp, "%+.16f ", u[t.database.ptr[0][i]+imove]); // my variable field
					avg_stat += u[t.database.ptr[0][i] + imove];
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			avg_stat /= t.database.maxelm;
			printf("average statistics=%e\n",avg_stat);
			/*
			if (bextendedprint) {
				// Pressure
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}
			*/
			fprintf(fp, "\n");

			/*
			// PAM
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// PAM
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				}
				else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Rho
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mu
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // Турбулентная динамическая вязкость
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", 1.0*f[t.ptr[1][i]].icolor_different_fluid_domain[t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// MUT
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == ZEROEQMOD) || (f[t.ptr[1][i]].iflowregime == SMAGORINSKY)) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime == ZEROEQMOD) || (f[0].iflowregime == SMAGORINSKY)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // Завихрённость - модуль ротора скорости.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Curl
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i + maxelm]); // GRADZVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");
			*/
		}
		/*
		doublereal *Tx = nullptr;
		doublereal *Ty = nullptr;
		doublereal *Tz = nullptr;
		Tx = new doublereal[t.maxelm + t.maxbound];
		Ty = new doublereal[t.maxelm + t.maxbound];
		Tz = new doublereal[t.maxelm + t.maxbound];

		// инициализация нулём.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// нахождение градиентов.
		for (i = 0; i<t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, false,
				t.sosedb, Tx, Ty, Tz);
		}

		for (i = 0; i<t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, true,
				t.sosedb, Tx, Ty, Tz);
		}


		// Сохранение в файл.

		// Heat Flux X
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", -t.prop[LAM][i] * Tx[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Heat Flux Y
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", -t.prop[LAM][i] * Ty[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Ty[i + maxelm]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Heat Flux Z
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", -t.prop[LAM][i] * Tz[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tz[i + maxelm]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Mag Heat Flux
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i])));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
				fprintf(fp, "%+.16f ", sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm])));

				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");


		// Освобождение оперативной памяти.
		if (Tx != nullptr) {
			delete[] Tx;
		}
		if (Ty != nullptr) {
			delete[] Ty;
		}
		if (Tz != nullptr) {
			delete[] Tz;
		}
		*/


		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {

#if doubleintprecision == 1
					// Визуализация твёрдого тела и только как и раньше.
					fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
					// Визуализация твёрдого тела и только как и раньше.
					fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;
						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {

#if doubleintprecision == 1
							fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

								}
					}
					else {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и жидкости.
						fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и жидкости.
						fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						}
				}
			}
		}
		else {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1 = fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
			if (fp1 == nullptr) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}
	//getchar();
	system("PAUSE");
	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

} // for amg


// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2_rev(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, BLOCK* &b, integer lb)
{
	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // часть 1 или 3
	errno_t err=0;
#ifdef MINGW_COMPILLER
	err = 0;
	fp = fopen64("ALICEFLOW0_07_temp.PLT", "w");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
#endif
	
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.


	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	// w -write, b - binary.
	if ((err) != 0) {
		printf("Create File temp Error in exporttecplotxy360T_3D_part2_rev in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1 = fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1 == NULL) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
		if (bOk)
		{


			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // считаем что температура всегда активна
				}
			}

			if (ivarexport == 1) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// запись имён переменных
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // печатается только поле температур
				fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, heat_flux_x, heat_flux_y, heat_flux_z, mag_heat_flux\n");  // печатается только поле температур


#if doubleintprecision == 1																													// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
			    }
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#else																											// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#endif

				

				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			else if (ivarexport == 3) {
				// запись заголовка
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer ncell_shadow = ncell;
				if ((ionly_solid_visible == 1) && (flow_interior > 0))
				{
					ncell_shadow = 0;
					integer ncell_shadow_old = 0;
					for (i = 0; i < t.database.ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i]-1;
						inode2 = t.database.nvtxcell[1][i]-1;
						inode3 = t.database.nvtxcell[2][i]-1;
						inode4 = t.database.nvtxcell[3][i]-1;
						inode5 = t.database.nvtxcell[4][i]-1;
						inode6 = t.database.nvtxcell[5][i]-1;
						inode7 = t.database.nvtxcell[6][i]-1;
						inode8 = t.database.nvtxcell[7][i]-1;
						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
						&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {
						    ncell_shadow++;
						}
						/*
						// Верен и работоспособен вышеприведённый код.
						integer ib;
						TOCHKA p;
						p.x = 0.5*(t.pa[t.nvtx[1][inode1] - 1].x + t.pa[t.nvtx[0][inode1] - 1].x);
						p.y = 0.5*(t.pa[t.nvtx[2][inode1] - 1].y + t.pa[t.nvtx[0][inode1] - 1].y);
						p.z = 0.5*(t.pa[t.nvtx[4][inode1] - 1].z + t.pa[t.nvtx[0][inode1] - 1].z);
						if (!in_model_flow(p, ib, b, lb)) {
							p.x = 0.5*(t.pa[t.nvtx[1][inode2] - 1].x + t.pa[t.nvtx[0][inode2] - 1].x);
							p.y = 0.5*(t.pa[t.nvtx[2][inode2] - 1].y + t.pa[t.nvtx[0][inode2] - 1].y);
							p.z = 0.5*(t.pa[t.nvtx[4][inode2] - 1].z + t.pa[t.nvtx[0][inode2] - 1].z);
							if (!in_model_flow(p, ib, b, lb)) {
								p.x = 0.5*(t.pa[t.nvtx[1][inode3] - 1].x + t.pa[t.nvtx[0][inode3] - 1].x);
								p.y = 0.5*(t.pa[t.nvtx[2][inode3] - 1].y + t.pa[t.nvtx[0][inode3] - 1].y);
								p.z = 0.5*(t.pa[t.nvtx[4][inode3] - 1].z + t.pa[t.nvtx[0][inode3] - 1].z);
								if (!in_model_flow(p, ib, b, lb)) {
									p.x = 0.5*(t.pa[t.nvtx[1][inode4] - 1].x + t.pa[t.nvtx[0][inode4] - 1].x);
									p.y = 0.5*(t.pa[t.nvtx[2][inode4] - 1].y + t.pa[t.nvtx[0][inode4] - 1].y);
									p.z = 0.5*(t.pa[t.nvtx[4][inode4] - 1].z + t.pa[t.nvtx[0][inode4] - 1].z);
									if (!in_model_flow(p, ib, b, lb)) {
										p.x = 0.5*(t.pa[t.nvtx[1][inode5] - 1].x + t.pa[t.nvtx[0][inode5] - 1].x);
										p.y = 0.5*(t.pa[t.nvtx[2][inode5] - 1].y + t.pa[t.nvtx[0][inode5] - 1].y);
										p.z = 0.5*(t.pa[t.nvtx[4][inode5] - 1].z + t.pa[t.nvtx[0][inode5] - 1].z);
										if (!in_model_flow(p, ib, b, lb)) {
											p.x = 0.5*(t.pa[t.nvtx[1][inode6] - 1].x + t.pa[t.nvtx[0][inode6] - 1].x);
											p.y = 0.5*(t.pa[t.nvtx[2][inode6] - 1].y + t.pa[t.nvtx[0][inode6] - 1].y);
											p.z = 0.5*(t.pa[t.nvtx[4][inode6] - 1].z + t.pa[t.nvtx[0][inode6] - 1].z);
											if (!in_model_flow(p, ib, b, lb)) {
												p.x = 0.5*(t.pa[t.nvtx[1][inode7] - 1].x + t.pa[t.nvtx[0][inode7] - 1].x);
												p.y = 0.5*(t.pa[t.nvtx[2][inode7] - 1].y + t.pa[t.nvtx[0][inode7] - 1].y);
												p.z = 0.5*(t.pa[t.nvtx[4][inode7] - 1].z + t.pa[t.nvtx[0][inode7] - 1].z);
												if (!in_model_flow(p, ib, b, lb)) {
													p.x = 0.5*(t.pa[t.nvtx[1][inode8] - 1].x + t.pa[t.nvtx[0][inode8] - 1].x);
													p.y = 0.5*(t.pa[t.nvtx[2][inode8] - 1].y + t.pa[t.nvtx[0][inode8] - 1].y);
													p.z = 0.5*(t.pa[t.nvtx[4][inode8] - 1].z + t.pa[t.nvtx[0][inode8] - 1].z);
													if (!in_model_flow(p, ib, b, lb)) {
														ncell_shadow++;
													}
												}
											}
										}
									}
								}
							}
						}
						*/
					}

#if doubleintprecision == 1
					//printf("%lld %lld\n", ncell_shadow, ncell_shadow_old);
#else
					//printf("%d %d\n", ncell_shadow, ncell_shadow_old);
#endif

					
					//getchar();
					if (ncell_shadow == 0) {
						// У нас совсем нету твёрдого тела, поэтому мы показываем только  жидкость.
						ncell_shadow = ncell;
						ionly_solid_visible = 0;
					}
				}
				// Полный набор искомых величин и теплопередача и гидродинамика:
				if (bextendedprint) {
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
				}
				else {
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
				}

#if doubleintprecision == 1
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
			    }
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// запись информации о зонах
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif

				
				if (bvery_big_memory) {
					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			if (!bvery_big_memory) {
				fclose(fp1); // закрытие файла
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части

		// Запись поля температур производится всегда.


		// запись температуры
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.potent[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = maxelm; i<maxelm + t.maxbound; i++) {
				fprintf(fp, "%+.16f ", t.potent[i]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}


		fprintf(fp, "\n");

		// Lam
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.prop[LAM][i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Запись гидродинамических величин если необходимо:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VX][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VY][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VX][i + maxelm] * f[idfluid].potent[VX][i + maxelm];
					doublereal svy = f[idfluid].potent[VY][i + maxelm] * f[idfluid].potent[VY][i + maxelm];
					doublereal svz = f[idfluid].potent[VZ][i + maxelm] * f[idfluid].potent[VZ][i + maxelm];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}


			fprintf(fp, "\n");

			// Pressure
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Pressure
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// PAM
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// PAM
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				}
				else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Rho
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mu
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // Турбулентная динамическая вязкость
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// MUT
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == ZEROEQMOD) || 
						(f[t.ptr[1][i]].iflowregime == SMAGORINSKY) ||
						(f[t.ptr[1][i]].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime == ZEROEQMOD) || 
						(f[0].iflowregime == SMAGORINSKY)||
						(f[0].iflowregime == RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == RANS_MENTER_SST) ||
						(f[0].iflowregime == RANS_STANDART_K_EPS)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // Завихрённость - модуль ротора скорости.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Curl
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVX
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADXVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADYVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// GRADZVZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i + maxelm]); // GRADZVZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

		}

		doublereal *Tx = nullptr;
		doublereal *Ty = nullptr;
		doublereal *Tz = nullptr;
		Tx = new doublereal[t.maxelm + t.maxbound];
		Ty = new doublereal[t.maxelm + t.maxbound];
		Tz = new doublereal[t.maxelm + t.maxbound];

		// инициализация нулём.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// нахождение градиентов.
		for (i = 0; i<t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, false,
				t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
		}

		for (i = 0; i<t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.sosedi, t.maxelm, true,
				t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
		}


		// Сохранение в файл.

		// Heat Flux X
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", -t.prop[LAM][i] * Tx[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Heat Flux Y
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", -t.prop[LAM][i] * Ty[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Ty[i + maxelm]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Heat Flux Z
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", -t.prop[LAM][i] * Tz[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tz[i + maxelm]);
				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");

		// Mag Heat Flux
		for (i = 0; i<maxelm; i++) {
			fprintf(fp, "%+.16f ", sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i])));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i = 0; i<t.maxbound; i++) {
				//fprintf(fp, "%+.16f ", -t.prop_b[LAM][i] * Tx[i + maxelm]);
				fprintf(fp, "%+.16f ", sqrt((-t.prop_b[LAM][i] * Tx[i + maxelm])*(-t.prop_b[LAM][i] * Tx[i + maxelm]) + (-t.prop_b[LAM][i] * Ty[i + maxelm])*(-t.prop_b[LAM][i] * Ty[i + maxelm]) + (-t.prop_b[LAM][i] * Tz[i + maxelm])*(-t.prop_b[LAM][i] * Tz[i + maxelm])));

				if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n");


		// Освобождение оперативной памяти.
		if (Tx != nullptr) {
			delete[] Tx;
		}
		if (Ty != nullptr) {
			delete[] Ty;
		}
		if (Tz != nullptr) {
			delete[] Tz;
		}



		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if ((ionly_solid_visible == 1) && (flow_interior>0))
				{
					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i]-1;
					inode2 = t.database.nvtxcell[1][i]-1;
					inode3 = t.database.nvtxcell[2][i]-1;
					inode4 = t.database.nvtxcell[3][i]-1;
					inode5 = t.database.nvtxcell[4][i]-1;
					inode6 = t.database.nvtxcell[5][i]-1;
					inode7 = t.database.nvtxcell[6][i]-1;
					inode8 = t.database.nvtxcell[7][i]-1;
					// визуализация только твёрдого тела.
					// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
					if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
						(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {

#if doubleintprecision == 1
						fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						}
				}
				else {


#if doubleintprecision == 1
					// Визуализация твёрдого тела и жидкости.
					fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
					// Визуализация твёрдого тела и жидкости.
					fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

					}
			}
		}
		else {
#ifdef MINGW_COMPILLER
			err = 0;
			fp1 = fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
			if (fp1 == nullptr) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != nullptr) {
					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// Файл не был открыт.
		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}//30 april 2016


// Это аналог VARIATION Plot из ANSYS Icepak.
// С 16 января 2018 года мы передаём из интерфейса АЛИСМеш_v0_42 информацию о
// линии на которой хотим отобразить температуру (передаётся опорная точка через
// которую линия точно проходит и направление линии, которое должно совпадать с направлением одной
// из осей декартовой прямоугольной системы координат.
void xyplot( FLOW* &fglobal, integer flow_interior, TEMPER &t) {
	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	err = 0;
	fp = fopen64("xyplot.txt", "w");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "xyplot.txt", "w");
#endif
	

	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else {

		TOCHKA p;
		doublereal epsilon=1e30;
		doublereal dist;
		doublereal x=0.0e-3, y=0.0e-3, z=0.0e-3; // точка через которую проходит линия
		x = Tochka_position_X0_for_XY_Plot;
		y = Tochka_position_Y0_for_XY_Plot;
		z = Tochka_position_Z0_for_XY_Plot;

		integer ifi=0, iPf=0;
		integer iplane=XZ; // плоскость перпендикулярная линии.
		switch (idirectional_for_XY_Plot) {
			case 0 : // X
			    iplane = YZ; // плоскость перпендикулярная линии.
			    break;
			case 1 : // Y
			    iplane = XZ; // плоскость перпендикулярная линии.
		    	break;
			case 2 : // Z
			    iplane = XY; // плоскость перпендикулярная линии.
			    break;
			default : // X
				iplane = YZ; // плоскость перпендикулярная линии.
				break;
		}


		for (integer i=0; i<flow_interior; i++) {
			for (integer iP=0; iP<fglobal[i].maxelm; iP++) {
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p,100); // вычисление координат центра КО.
				dist=sqrt(fabs(x-p.x)*fabs(x-p.x)+fabs(y-p.y)*fabs(y-p.y)+fabs(z-p.z)*fabs(z-p.z));
				if (dist<epsilon) {
					epsilon=dist;
					ifi=i;
					iPf=iP;
				}
			}
		}
		// перемотка в начало 
		switch (iplane) {
		case XY: while (fglobal[ifi].sosedi[BSIDE][iPf].iNODE1<fglobal[ifi].maxelm) iPf = fglobal[ifi].sosedi[BSIDE][iPf].iNODE1; break;
		case XZ: while (fglobal[ifi].sosedi[SSIDE][iPf].iNODE1<fglobal[ifi].maxelm) iPf = fglobal[ifi].sosedi[SSIDE][iPf].iNODE1; break;
		case YZ: while (fglobal[ifi].sosedi[WSIDE][iPf].iNODE1<fglobal[ifi].maxelm) iPf = fglobal[ifi].sosedi[WSIDE][iPf].iNODE1; break;
		}

		integer G;
		switch (iplane) {
		  case XY : G=TSIDE;  break;
		  case XZ : G=NSIDE;  break;
		  case YZ : G=ESIDE;  break;
		}
		
		
		fprintf(fp, "position,\tVx,\tVy,\tVz,\tPam,\tPress,\tFbuf,\tGRADPRESS,\tGRADPAM,\tmassfluxingran\n");
		doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	    volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p,100); // вычисление координат центра КО.
		switch (iplane) {
		  case XY : fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			  p.z - 0.5*dz, fglobal[ifi].potent[VX][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[VY][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[VZ][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[PAM][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[PRESS][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[FBUF][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[GRADZPRESS][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[GRADZPAM][fglobal[ifi].sosedi[BSIDE][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
			              break;
		  case XZ :  fprintf(fp,"%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			  p.y - 0.5*dy, fglobal[ifi].potent[VX][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[VY][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[VZ][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[PAM][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[PRESS][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[FBUF][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[GRADYPRESS][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[GRADYPAM][fglobal[ifi].sosedi[SSIDE][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
			              break;
		  case YZ :  fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			  p.x - 0.5*dx, fglobal[ifi].potent[VX][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[VY][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[VZ][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[PAM][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[PRESS][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[FBUF][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[GRADXPRESS][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
			  fglobal[ifi].potent[GRADXPAM][fglobal[ifi].sosedi[WSIDE][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
			              break;
		}
		switch (iplane) {
		  case XY : while (iPf<fglobal[ifi].maxelm) {
			        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p,100); 
			        fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n", 
						  p.z, fglobal[ifi].potent[VX][iPf],
					      fglobal[ifi].potent[VY][iPf],
						  fglobal[ifi].potent[VZ][iPf],
						  fglobal[ifi].potent[PAM][iPf],
						  fglobal[ifi].potent[PRESS][iPf],
						  fglobal[ifi].potent[FBUF][iPf],
						  fglobal[ifi].potent[GRADZPRESS][iPf],
						  fglobal[ifi].potent[GRADZPAM][iPf],
						  fglobal[ifi].mf[iPf][G]);
					if (fglobal[ifi].sosedi[TSIDE][iPf].iNODE1 >= fglobal[ifi].maxelm) {
						  volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
							  p.z + 0.5*dz, fglobal[ifi].potent[VX][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[VY][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[VZ][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[PAM][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[PRESS][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[FBUF][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[GRADZPRESS][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[GRADZPAM][fglobal[ifi].sosedi[TSIDE][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
					}
					iPf = fglobal[ifi].sosedi[TSIDE][iPf].iNODE1;
					} break;
		  case XZ : while (iPf<fglobal[ifi].maxelm) {
			        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p,100); 
			        fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f  %+.16f %+.16f %+.16f\n", 
						  p.y, fglobal[ifi].potent[VX][iPf],
					      fglobal[ifi].potent[VY][iPf],
						  fglobal[ifi].potent[VZ][iPf],
						  fglobal[ifi].potent[PAM][iPf],
						  fglobal[ifi].potent[PRESS][iPf],
						  fglobal[ifi].potent[FBUF][iPf],
						  fglobal[ifi].potent[GRADYPRESS][iPf],
						  fglobal[ifi].potent[GRADYPAM][iPf],
						  fglobal[ifi].mf[iPf][G]);
					if (fglobal[ifi].sosedi[NSIDE][iPf].iNODE1 >= fglobal[ifi].maxelm) {
						  volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
							  p.y + 0.5*dy, fglobal[ifi].potent[VX][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[VY][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[VZ][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[PAM][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[PRESS][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[FBUF][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[GRADYPRESS][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
							  fglobal[ifi].potent[GRADYPAM][fglobal[ifi].sosedi[NSIDE][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
					}
					/*

					#if doubleintprecision == 1
						// Узнаём последовательность узлов для отладки.
						printf("iPf=%lld\n",iPf);
						if (fglobal[ifi].sosedi[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
						    printf("iPffinish=%lld\n",fglobal[ifi].sosedi[NSIDE][iPf].iNODE1);
							getchar();
						}
					#else
						// Узнаём последовательность узлов для отладки.
						printf("iPf=%d\n",iPf);
						if (fglobal[ifi].sosedi[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
							printf("iPffinish=%d\n",fglobal[ifi].sosedi[NSIDE][iPf].iNODE1);
							getchar();
						}
					#endif

					
					*/
					iPf = fglobal[ifi].sosedi[NSIDE][iPf].iNODE1;
					} break;
		  case YZ : while (iPf<fglobal[ifi].maxelm) {
			        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p,100); 
			        fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n", 
						  p.y, fglobal[ifi].potent[VX][iPf],
					      fglobal[ifi].potent[VY][iPf],
						  fglobal[ifi].potent[VZ][iPf],
						  fglobal[ifi].potent[PAM][iPf],
						  fglobal[ifi].potent[PRESS][iPf],
						  fglobal[ifi].potent[FBUF][iPf],
						  fglobal[ifi].potent[GRADXPRESS][iPf],
						  fglobal[ifi].potent[GRADXPAM][iPf],
						  fglobal[ifi].mf[iPf][G]);
					if (fglobal[ifi].sosedi[ESIDE][iPf].iNODE1 >= fglobal[ifi].maxelm) {
						  volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
							  p.x + 0.5*dx, fglobal[ifi].potent[VX][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[VY][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[VZ][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[PAM][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[PRESS][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[FBUF][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[GRADXPRESS][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
							  fglobal[ifi].potent[GRADXPAM][fglobal[ifi].sosedi[ESIDE][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
					}
					iPf = fglobal[ifi].sosedi[ESIDE][iPf].iNODE1;
					} break;
		}
		fclose(fp);
	}
}

// построение графика температуры вдоль линии 
void xyplot_temp(TEMPER &t, doublereal* tempfiltr) {

	// tempfiltr - передаваемая однократно фильтрованная температура.
	// имя создаваемого файла xyplotT.txt.

	FILE *fp=NULL;
	errno_t err=0;
#ifdef MINGW_COMPILLER
	fp = fopen64("xyplotT.txt", "w");
	if (fp == NULL) err = 1;
#else
	err = fopen_s(&fp, "xyplotT.txt", "w");
#endif
	

	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else if (t.maxelm > 0) {
		if (fp != NULL)
		{

			// внимание ! требуется указать точку через которую будет проходить линия и плоскость которой данная линия будет перпендикулярна.
			TOCHKA p;
			p.x = 0;
			p.y = 0;
			p.z = 0;
			doublereal epsilon = 1e30;
			doublereal dist;
			//doublereal x = 2.0245e-3, y = 0.2268e-3, z = 0.0125e-3; // точка через которую проходит линия
			doublereal x = 3e-3, y = -0.00005e-3, z = 0.57425e-3;
			x = Tochka_position_X0_for_XY_Plot;
			y = Tochka_position_Y0_for_XY_Plot;
			z = Tochka_position_Z0_for_XY_Plot;
			integer iPf = 0;
			integer iplane = YZ; // плоскость перпендикулярная линии.
			switch (idirectional_for_XY_Plot) {
			case 0: // YZ 
				iplane = YZ; // плоскость перпендикулярная линии.
				break;
			case 1: //XZ
				iplane = XZ; // плоскость перпендикулярная линии.
				break;
			case 2: // XY
				iplane = XY; // плоскость перпендикулярная линии.
				break;
			default:
				iplane = YZ; // плоскость перпендикулярная линии.
				break;
			}

			for (integer iP = 0; iP < t.maxelm; iP++) {
				center_cord3D(iP, t.nvtx, t.pa, p,100); // вычисление координат центра КО.
				dist = sqrt(fabs(x - p.x)*fabs(x - p.x) + fabs(y - p.y)*fabs(y - p.y) + fabs(z - p.z)*fabs(z - p.z));
				if (dist < epsilon) {
					epsilon = dist;
					iPf = iP;
				}
			}

			// перемотка в начало 
			if (t.maxelm > 0) {
				switch (iplane) {
				case XY: while (t.sosedi[BSIDE][iPf].iNODE1 < t.maxelm) iPf = t.sosedi[BSIDE][iPf].iNODE1; break;
				case XZ: while (t.sosedi[SSIDE][iPf].iNODE1 < t.maxelm) iPf = t.sosedi[SSIDE][iPf].iNODE1; break;
				case YZ: while (t.sosedi[WSIDE][iPf].iNODE1 < t.maxelm) iPf = t.sosedi[WSIDE][iPf].iNODE1; break;
				}
			}

			integer G;
			switch (iplane) {
			case XY: G = TSIDE;  break;
			case XZ: G = NSIDE;  break;
			case YZ: G = ESIDE;  break;
			}


			fprintf(fp, "position,\ttemperature,\ttemperature_avg\n");
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
			if (t.maxelm > 0) {
				volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
				center_cord3D(iPf, t.nvtx, t.pa, p, 100); // вычисление координат центра КО.
			}
			switch (iplane) {
			case XY: fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.z - 0.5*dz, t.potent[t.sosedi[BSIDE][iPf].iNODE1],
				tempfiltr[t.sosedi[BSIDE][iPf].iNODE1]);
				break;
			case XZ:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.y - 0.5*dy, t.potent[t.sosedi[SSIDE][iPf].iNODE1],
				tempfiltr[t.sosedi[SSIDE][iPf].iNODE1]);
				break;
			case YZ:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.x - 0.5*dx, t.potent[t.sosedi[WSIDE][iPf].iNODE1],
				tempfiltr[t.sosedi[WSIDE][iPf].iNODE1]);
				break;
			}
			switch (iplane) {
			case XY: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p,100);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.z, t.potent[iPf],
					tempfiltr[iPf]);
				if (t.sosedi[TSIDE][iPf].iNODE1 >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.z + 0.5*dz, t.potent[t.sosedi[TSIDE][iPf].iNODE1],
						tempfiltr[t.sosedi[TSIDE][iPf].iNODE1]);
				}
				iPf = t.sosedi[TSIDE][iPf].iNODE1;
			} break;
			case XZ: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p,100);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.y, t.potent[iPf],
					tempfiltr[iPf]);
				if (t.sosedi[NSIDE][iPf].iNODE1 >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.y + 0.5*dy, t.potent[t.sosedi[NSIDE][iPf].iNODE1],
						tempfiltr[t.sosedi[NSIDE][iPf].iNODE1]);
				}
				/*

				#if doubleintprecision == 1
					// Узнаём последовательность узлов для отладки.
					printf("iPf=%lld\n",iPf);
					if (fglobal[ifi].sosedi[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
						printf("iPffinish=%lld\n",fglobal[ifi].sosedi[NSIDE][iPf].iNODE1);
						getchar();
					}
				#else
					// Узнаём последовательность узлов для отладки.
					printf("iPf=%d\n",iPf);
					if (fglobal[ifi].sosedi[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
						printf("iPffinish=%d\n",fglobal[ifi].sosedi[NSIDE][iPf].iNODE1);
						getchar();
					}
				#endif

				
				*/
				iPf = t.sosedi[NSIDE][iPf].iNODE1;
			} break;
			case YZ: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p,100);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.x, t.potent[iPf],
					tempfiltr[iPf]);
				if (t.sosedi[ESIDE][iPf].iNODE1 >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.x + 0.5*dx, t.potent[t.sosedi[ESIDE][iPf].iNODE1],
						tempfiltr[t.sosedi[ESIDE][iPf].iNODE1]);
				}
				iPf = t.sosedi[ESIDE][iPf].iNODE1;
			} break;
			}
			fclose(fp);
		}
	}
} // xyplot_temp

void animationtecplot360T_3D_part2(integer maxelm, integer ncell,
	FLOW* &f, TEMPER &t, integer flow_interior_count, char* title,
	bool btitle, integer iVar, BLOCK* &b, integer &lb)
{
    // title - заголовок зоны.
	// btitle - печатать ли общий заголовок

	FILE *fp=NULL;
    FILE *fp1=NULL; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	
	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	//
#ifdef MINGW_COMPILLER
	err = 0;
	if (btitle) {
		// мы стираем предыдущие кадры и переходим к новой анимации:
		switch (iVar) {
		case TEMP:  fp = fopen64("ALICEFLOW0_07_animation_temp.PLT", "w");
			break;
		case SPEED: fp = fopen64("ALICEFLOW0_07_animation_speed.PLT", "w");
			break;
		case PRESS: fp = fopen64("ALICEFLOW0_07_animation_press.PLT", "w");
			break;
		case PAM:  fp = fopen64("ALICEFLOW0_07_animation_pam.PLT", "w");
			break;
		}
	}
	else {
		// мы добавляем следующие анимационные кадры.
		switch (iVar) {
		case TEMP: fp = fopen64("ALICEFLOW0_07_animation_temp.PLT", "a");
			break;
		case SPEED: fp  = fopen64("ALICEFLOW0_07_animation_speed.PLT", "a");
			break;
		case PRESS: fp = fopen64("ALICEFLOW0_07_animation_press.PLT", "a");
			break;
		case PAM:  fp = fopen64("ALICEFLOW0_07_animation_pam.PLT", "a");
			break;
		}
	}
	if (fp == NULL) err = 1;
#else
	if (btitle) {
		// мы стираем предыдущие кадры и переходим к новой анимации:
		switch (iVar) {
		case TEMP: err = fopen_s(&fp, "ALICEFLOW0_07_animation_temp.PLT", "w");
			break;
		case SPEED: err = fopen_s(&fp, "ALICEFLOW0_07_animation_speed.PLT", "w");
			break;
		case PRESS: err = fopen_s(&fp, "ALICEFLOW0_07_animation_press.PLT", "w");
			break;
		case PAM:  err = fopen_s(&fp, "ALICEFLOW0_07_animation_pam.PLT", "w");
			break;
		}
	}
	else {
		// мы добавляем следующие анимационные кадры.
		switch (iVar) {
		case TEMP: err = fopen_s(&fp, "ALICEFLOW0_07_animation_temp.PLT", "a");
			break;
		case SPEED: err = fopen_s(&fp, "ALICEFLOW0_07_animation_speed.PLT", "a");
			break;
		case PRESS: err = fopen_s(&fp, "ALICEFLOW0_07_animation_press.PLT", "a");
			break;
		case PAM:  err = fopen_s(&fp, "ALICEFLOW0_07_animation_pam.PLT", "a");
			break;
		}
	}
#endif
	
	
	if (err != 0) {
		printf("Create File temp Error in animationtecplot360T_3D_part2 in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {
        
       // int c; // читаемый символ
		integer ivarexport=1; // по умолчанию только поле температур:
		integer i=0; // счётчик цикла
		/*
#ifdef MINGW_COMPILLER
		err = 0;
		fp1 = fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
		if (fp1 == NULL) err = 1;
#else
		err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif
        if ((err) != 0) {
		    printf("Open File temp part1 Error\n");
		    //getchar();
			system("pause");

	    }
	    else 
		*/
		{
			// копирование первой части в итоговый файл
			// Особенность : иногда необходимо изменить вторую строку в файле:
			if (flow_interior_count>0) {
				// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
				for (i=0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport=3; // считаем что температура всегда активна
				}
			}
            
			if (ivarexport==1) {

				// запись заголовка
		        fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

		         // запись имён переменных
				switch (iVar) {
	               case TEMP : fprintf(fp, "VARIABLES = x, y, z, Temp\n");  // печатается только поле температур
		                       break;
	               case SPEED :  fprintf(fp, "VARIABLES = x, y, z, Speed\n"); // печатается модуль скорости
		                       break;
	               case PRESS :  fprintf(fp, "VARIABLES = x, y, z, Press\n"); // печатается давление
		                       break;
	               case PAM :   fprintf(fp, "VARIABLES = x, y, z, PAM\n"); // печатается поправка давления
		                       break;
	            }
		         
#if doubleintprecision == 1
				// запись информации о зонах
				fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
				// запись информации о зонах
				fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif
				
		
				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (1 && 2 == steady_or_unsteady_global_determinant) {
						// quolity
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
							volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", dmax / dmin);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						// log10(quolity)
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
							volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", log10(dmax / dmin));
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				//while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

			}
			else if (ivarexport==3) {
				
				if (btitle) {

				    // запись заголовка
		            fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

     				// запись имён переменных
				    switch (iVar) {
	                   case TEMP : fprintf(fp, "\nVARIABLES = x, y, z, Temp\n");  // печатается только поле температур
		                           break;
	                   case SPEED :  fprintf(fp, "\nVARIABLES = x, y, z, Speed\n"); // печатается модуль скорости
		                           break;
	                   case PRESS :  fprintf(fp, "\nVARIABLES = x, y, z, Press\n"); // печатается давление
		                           break;
	                   case PAM :   fprintf(fp, "\nVARIABLES = x, y, z, PAM\n"); // печатается поправка давления
		                           break;
	                }
				}
					
				// запись информации о зонах
#if doubleintprecision == 1
				fprintf(fp, "ZONE T=\"%s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", title, maxelm, ncell);
#else
				fprintf(fp, "ZONE T=\"%s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", title, maxelm, ncell);
#endif
                
				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (1 && 2 == steady_or_unsteady_global_determinant) {
						// quolity
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
							volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", dmax / dmin);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						// log10(quolity)
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
							volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", log10(dmax / dmin));
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				//while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

						
			}
            //fclose(fp1); // закрытие файла
            //printf("export tecplot part1 is successfully reading and written...OK.\n");
		}

		// запись второй части
        
        // Запись поля температур производится всегда.

        if (iVar==TEMP) {
			// запись температуры
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", t.potent[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

    		fprintf(fp, "\n");
		}
		

		/*
		// Lam
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.prop[LAM][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");
		*/

		// Запись гидродинамических величин если необходимо:
		if (ivarexport==3) {

			if (iVar==SPEED) {
				// Speed
                for (i=0;i<maxelm; i++) {
				    if (t.ptr[1][i]>-1) {
					    doublereal svx=f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					    doublereal svy=f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					    doublereal svz=f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					    fprintf(fp, "%+.16f ",sqrt(svx+svy+svz)); // f[t.ptr[1][i]].potent[FBUF][i]
				    } else fprintf(fp, "%+.16f ", 0.0);
                    if (i%10==0) fprintf(fp, "\n");
		        }

		        fprintf(fp, "\n");
			}

			if (iVar==PRESS) {
				// Pressure
                for (i=0;i<maxelm; i++) {
                    if (t.ptr[1][i]>-1) {
                       fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
			     	} else fprintf(fp, "%+.16f ", 0.0);
                    if (i%10==0) fprintf(fp, "\n");
		        }

		        fprintf(fp, "\n");
			}
			

			if (iVar==PAM) {
				// PAM
                for (i=0;i<maxelm; i++) {
                    if (t.ptr[1][i]>-1) {
                       fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				    } else fprintf(fp, "%+.16f ", 0.0);
                    if (i%10==0) fprintf(fp, "\n");
		        }

		        fprintf(fp, "\n");
			}
			
			/*
			// VX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/
			/*
			// VY
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/
			/*
			// VZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/
			/*
            // Rho
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				} else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/

			/*
			// Mu
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				} else fprintf(fp, "%+.16f ", 0.0); 
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/

			/*
			// Mut // Турбулентная динамическая вязкость
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/

			/*
			// для отладки график распределения нумерации контрольных объёмов.
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", doublereal(i));
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/

		}

		/*
#ifdef MINGW_COMPILLER
		err = 0;
		fp1 = fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
		if (fp1 == nullptr) err = 1;
#else
		err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
        if ((err) != 0) {
		    printf("Open File temp part3 Error\n");
		    //getchar();
			system("pause");

	    }
	    else {
			if (fp1 != nullptr) {
				// копирование третьей части в итоговый файл
				while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				fclose(fp1); // закрытие файла
				//printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}
		*/

		//if (bvery_big_memory) 
		{
			// запись информации о разностной сетке
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = t.database.nvtxcell[0][i] - 1;
					inode2 = t.database.nvtxcell[1][i] - 1;
					inode3 = t.database.nvtxcell[2][i] - 1;
					inode4 = t.database.nvtxcell[3][i] - 1;
					inode5 = t.database.nvtxcell[4][i] - 1;
					inode6 = t.database.nvtxcell[5][i] - 1;
					inode7 = t.database.nvtxcell[6][i] - 1;
					inode8 = t.database.nvtxcell[7][i] - 1;

					/*
					TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
					//TOCHKA pall;
					center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
					center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
					center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
					center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
					center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
					center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
					center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
					center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
					*/
					integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
					/*
					in_model_temp(p1, ib1, b, lb);
					in_model_temp(p2, ib2, b, lb);
					in_model_temp(p3, ib3, b, lb);
					in_model_temp(p4, ib4, b, lb);
					in_model_temp(p5, ib5, b, lb);
					in_model_temp(p6, ib6, b, lb);
					in_model_temp(p7, ib7, b, lb);
					in_model_temp(p8, ib8, b, lb);
					*/

					ib1 = t.whot_is_block[inode1];
					ib2 = t.whot_is_block[inode2];
					ib3 = t.whot_is_block[inode3];
					ib4 = t.whot_is_block[inode4];
					ib5 = t.whot_is_block[inode5];
					ib6 = t.whot_is_block[inode6];
					ib7 = t.whot_is_block[inode7];
					ib8 = t.whot_is_block[inode8];


					if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						integer inode2W = t.sosedi[WSIDE][inode1].iNODE1;
						integer inode3W = t.sosedi[WSIDE][inode4].iNODE1;
						integer inode6W = t.sosedi[WSIDE][inode5].iNODE1;
						integer inode7W = t.sosedi[WSIDE][inode8].iNODE1;


						integer inode5B = t.sosedi[BSIDE][inode1].iNODE1;
						integer inode6B = t.sosedi[BSIDE][inode2].iNODE1;
						integer inode7B = t.sosedi[BSIDE][inode3].iNODE1;
						integer inode8B = t.sosedi[BSIDE][inode4].iNODE1;



						integer inode3S = t.sosedi[SSIDE][inode2].iNODE1;
						integer inode4S = t.sosedi[SSIDE][inode1].iNODE1;
						integer inode7S = t.sosedi[SSIDE][inode6].iNODE1;
						integer inode8S = t.sosedi[SSIDE][inode5].iNODE1;


						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8, pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);
						pall.x = 0.125*(p1.x + p2.x + p3.x + p4.x + p5.x + p6.x + p7.x + p8.x);
						pall.y = 0.125*(p1.y + p2.y + p3.y + p4.y + p5.y + p6.y + p7.y + p8.y);
						pall.z = 0.125*(p1.z + p2.z + p3.z + p4.z + p5.z + p6.z + p7.z + p8.z);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);

						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

							}
						}
					}
					else {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						/*
						TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
						//TOCHKA pall;
						center_cord3D(inode1, t.nvtx, t.pa, p1, 100);
						center_cord3D(inode2, t.nvtx, t.pa, p2, 100);
						center_cord3D(inode3, t.nvtx, t.pa, p3, 100);
						center_cord3D(inode4, t.nvtx, t.pa, p4, 100);
						center_cord3D(inode5, t.nvtx, t.pa, p5, 100);
						center_cord3D(inode6, t.nvtx, t.pa, p6, 100);
						center_cord3D(inode7, t.nvtx, t.pa, p7, 100);
						center_cord3D(inode8, t.nvtx, t.pa, p8, 100);

						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						in_model_temp(p1, ib1, b, lb);
						in_model_temp(p2, ib2, b, lb);
						in_model_temp(p3, ib3, b, lb);
						in_model_temp(p4, ib4, b, lb);
						in_model_temp(p5, ib5, b, lb);
						in_model_temp(p6, ib6, b, lb);
						in_model_temp(p7, ib7, b, lb);
						in_model_temp(p8, ib8, b, lb);

						if ((ib1 != t.whot_is_block[inode1]) || (ib2 != t.whot_is_block[inode2]) || (ib3 != t.whot_is_block[inode3]) || (ib4 != t.whot_is_block[inode4]) ||
							(ib5 != t.whot_is_block[inode5]) || (ib6 != t.whot_is_block[inode6]) || (ib7 != t.whot_is_block[inode7]) || (ib8 != t.whot_is_block[inode8])) {
							printf("SIGNAL MESSAGE ZAMENA NEVERNA."); getchar();
						}
						*/
						integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
						ib1 = t.whot_is_block[inode1];
						ib2 = t.whot_is_block[inode2];
						ib3 = t.whot_is_block[inode3];
						ib4 = t.whot_is_block[inode4];
						ib5 = t.whot_is_block[inode5];
						ib6 = t.whot_is_block[inode6];
						ib7 = t.whot_is_block[inode7];
						ib8 = t.whot_is_block[inode8];

						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}
		}


		fclose(fp); // закрытие файла
		//printf("export tecplot is successfully written...OK.\n");
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

} // анимация модуля скорости

// анимация одной функции в программе tecplot 360
// из-за ограничений tecplot 360 по памяти имеет смысл загружать в неё
// не более одной трёхмерной функции за раз. Если требуется анимировать несколько функций, то
// для этих целей подходит создание нескольких файлов по одному на каждую функцию.
void animationtecplot360T_3D_part2all(integer maxelm, integer ncell,
	FLOW* &f, TEMPER &t, integer flow_interior_count, char* title, bool btitle,
	BLOCK* &b, integer &lb)
{
	// величины которые следует анимировать.
	bool bTEMP=true;
	bool bSpeed=true;
	bool bPressure=true;
	bool bPAM=true;

	// Непосредственно запись анимации в разные файлы.
	if (bTEMP) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, TEMP,b,lb);
	}

	if (bSpeed) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, SPEED,b,lb);
	}

	if (bPressure) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, PRESS, b, lb);
	}

	if (bPAM) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, PAM, b, lb);
	}
	
}

#endif

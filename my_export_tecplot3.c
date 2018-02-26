// Файл my_export_tecplot3.c передача результатов
// моделирования в программу tecplot360

#ifndef MY_EXPORT_TECPLOT3_C
#define MY_EXPORT_TECPLOT3_C 1

//#include "windows.h" // для функции WinExec
#include <string.h>

// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
void exporttecplotxy360_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal** potent, doublereal **rhie_chow)
{
	FILE *fp=NULL;
	errno_t err;
	err = fopen_s(&fp, "ALICEFLOW0_03.PLT", "w");
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

			if (fp != NULL) {
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
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, "ALICEFLOW0_03.PLT", "wb");
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
void exporttecplotxy360T_3D_part1and3(integer maxelm, integer maxbound, bool bextendedprint, integer ncell, 
									  integer** nvtx, integer** nvtxcell, TOCHKA* pa,
									  BOUND* sosedb, integer ivarexport, integer** ptr_out)
{

	if (bvery_big_memory) {

		database.maxelm = maxelm;
		database.ncell = ncell;

		// extended printeger не предусмотрено.

		// Если память уже выделялась ранее то её надо освободить.
		if (database.x != NULL) {
			free(database.x);
			database.x = NULL;
		}
		if (database.y != NULL) {
			free(database.y);
			database.y = NULL;
		}
		if (database.z != NULL) {
			free(database.z);
			database.z = NULL;
		}
		if (database.ptr != NULL) {
			for (integer i_92 = 0; i_92 < 2; i_92++) {
				delete[] database.ptr[i_92];
				database.ptr[i_92] = NULL;
			}
			delete[] database.ptr;
			database.ptr = NULL;
		}
		if (database.nvtxcell != NULL) {
			for (integer i_92 = 0; i_92 < 8; i_92++) {
				delete[] database.nvtxcell[i_92];
				database.nvtxcell[i_92] = NULL;
			}
			delete[] database.nvtxcell;
			database.nvtxcell = NULL;
		}

		database.x = (doublereal*)malloc(maxelm*sizeof(doublereal));
		database.y = (doublereal*)malloc(maxelm*sizeof(doublereal));
		database.z = (doublereal*)malloc(maxelm*sizeof(doublereal));
		database.ptr = new integer*[2];
		for (integer j = 0; j < 2; j++) {
			database.ptr[j] = new integer[maxelm];
		}
		for (integer j = 0; j < 2; j++) {
			for (integer i = 0; i < maxelm; i++) {
				database.ptr[j][i] = ptr_out[j][i];
			}
		}
		database.nvtxcell = new integer*[8];
		for (integer i = 0; i < 8; i++) {
			database.nvtxcell[i] = new integer[ncell];
		}

		// запись x
		for (integer i = 0; i < maxelm; i++) {
			database.x[i] = 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x);
		}
		// запись y
		for (integer i = 0; i < maxelm; i++) {
			database.y[i] = 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y);
		}

		// запись z
		for (integer i = 0; i < maxelm; i++) {
			database.z[i] = 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z);
		}

		// запись информации о разностной сетке
		for (integer i = 0; i < ncell; i++) {
			for (integer j = 0; j < 8; j++) {
				database.nvtxcell[j][i] = nvtxcell[j][i];
			}
		}

	}
	else {

		database.maxelm = 0;
		database.ncell = 0;
		//Запись в файл очень медленная.
		database.x = NULL;
		database.y = NULL;
		database.z = NULL;
		database.nvtxcell = NULL;
		database.ptr = NULL;


		// расширенная печать
		// При расширенной печати мы печатаем также и граничные узлы.
		// bextendedprint=true; расширенная печать. 

		// ivarexport == 1 печатается только поле температур,
		// ivarexport == 2 печатается только гидродинамика,
		// ivarexport == 3 печатается и поле температур и гидродинамика.

		FILE *fp;
		errno_t err;
		err = fopen_s(&fp, "ALICEFLOW0_06_temp_part1.txt", "w");
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

		if ((err = fopen_s(&fp, "ALICEFLOW0_06_temp_part3.txt", "w")) != 0) {
			printf("Create File temp part3 Error\n");
			//getchar();
			system("pause");
		}
		else {

			if (fp != NULL) {

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

	FILE *fp;
    FILE *fp1; // часть 1 или 3
	errno_t err;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");
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
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {
        
        char c; // читаемый символ
		integer ivarexport=1; // по умолчанию только поле температур:
		integer i=0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
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
					 fwrite(&number, sizeof(int), 1, fp);
					 char tit4[] = ", E=";//4
					 for (i = 0; i <= 3; i++) {
						 symbol = tit4[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					 fwrite(&ncell, sizeof(int), 1, fp);
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
					 fwrite(&number, sizeof(int), 1, fp);
					 char tit4[] = ", E=";//4
					 for (i = 0; i <= 3; i++) {
						 symbol = tit4[i];
						 fwrite(&symbol, sizeof(char), 1, fp);
					 }
					 fwrite(&ncell, sizeof(int), 1, fp);
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
					 for (i = 0; i < database.maxelm; i++) {
						 //fprintf(fp, "%+.16f ", database.x[i]);
						 fnumber = database.x[i];
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
					 for (i = 0; i < database.maxelm; i++) {
						 //fprintf(fp, "%+.16f ", database.y[i]);
						 fnumber = database.y[i];
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
					 for (i = 0; i < database.maxelm; i++) {
						 //fprintf(fp, "%+.16f ", database.z[i]);
						 fnumber = database.z[i];
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
				if (bextendedprint) {
					 fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z\n");
				}
				else {
                    fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z\n");
				}
					
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
					for (i = 0; i < database.maxelm; i++) {
						//fprintf(fp, "%+.16f ", database.x[i]);
						fnumber = database.x[i];
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
					for (i = 0; i < database.maxelm; i++) {
						//fprintf(fp, "%+.16f ", database.y[i]);
						fnumber = database.y[i];
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
					for (i = 0; i < database.maxelm; i++) {
						//fprintf(fp, "%+.16f ", database.z[i]);
						fnumber = database.z[i];
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
					if ((f[t.ptr[1][i]].iflowregime==ZEROEQMOD) || (f[t.ptr[1][i]].iflowregime==SMAGORINSKY)) {
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
					if ((f[0].iflowregime==ZEROEQMOD) || (f[0].iflowregime==SMAGORINSKY)) {
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

			doublereal *Tx=NULL;
            doublereal *Ty=NULL;
            doublereal *Tz=NULL;
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
					t.sosedb,  Tx, Ty, Tz);
			}

			for (i=0; i<t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, true, 
					t.sosedb,Tx, Ty, Tz);
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
			if (Tx != NULL) {
				delete[] Tx;
			}
			if (Ty != NULL) {
				delete[] Ty;
			}
			if (Tz != NULL) {
				delete[] Tz;
			}

		}

		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < database.ncell; i++) {
#if doubleintprecision == 1
				//fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
				//fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);

#else
				//fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
				//fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);

#endif
				integer number = database.nvtxcell[0][i];
				fwrite(&number, sizeof(int), 1, fp);
				char symbol = ' ';
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[1][i];
				fwrite(&number, sizeof(int), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[2][i];
				fwrite(&number, sizeof(int), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[3][i];
				fwrite(&number, sizeof(int), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[4][i];
				fwrite(&number, sizeof(int), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[5][i];
				fwrite(&number, sizeof(int), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[6][i];
				fwrite(&number, sizeof(int), 1, fp);
				fwrite(&symbol, sizeof(char), 1, fp);
				number = database.nvtxcell[7][i];
				fwrite(&number, sizeof(int), 1, fp);
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

				if (fp1 != NULL) {
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

void save_velocity_for_init(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count) {
	errno_t err_inicialization_data;
	FILE* fp_inicialization_data;
	err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "w");
	if (err_inicialization_data != 0) {
		// открытие неудачно или файл отсутствует.
		printf("Create File load.txt Error\n");
		//getchar();
		system("pause");
	}
	else {
#if doubleintprecision == 1
		fprintf_s(fp_inicialization_data, "%lld\n %lld\n", maxelm, ncell);
#else
		fprintf_s(fp_inicialization_data, "%d\n %d\n", maxelm, ncell);
#endif

		
		for (int i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf_s(fp_inicialization_data, "%e ", p.x);
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		for (int i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf_s(fp_inicialization_data, "%e ", p.y);
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");
		for (int i = 0; i < maxelm; i++) {
			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			fprintf_s(fp_inicialization_data, "%e ", p.z);
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");

		for (int i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VX][t.ptr[0][i]]);
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");

		for (int i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VY][t.ptr[0][i]]);
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");

		for (int i = 0; i < maxelm; i++) {
			if (t.ptr[1][i] > -1) {
				fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VZ][t.ptr[0][i]]);
			}
			else {
				fprintf_s(fp_inicialization_data, "%e ", 0.0);
			}
			if (i % 20 == 0) fprintf_s(fp_inicialization_data, "\n");
		}
		fprintf_s(fp_inicialization_data, "\n");

		for (int i = 0; i < ncell; i++) {
			integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
			inode1 = database.nvtxcell[0][i] - 1;
			inode2 = database.nvtxcell[1][i] - 1;
			inode3 = database.nvtxcell[2][i] - 1;
			inode4 = database.nvtxcell[3][i] - 1;
			inode5 = database.nvtxcell[4][i] - 1;
			inode6 = database.nvtxcell[5][i] - 1;
			inode7 = database.nvtxcell[6][i] - 1;
			inode8 = database.nvtxcell[7][i] - 1;
#if doubleintprecision == 1
			fprintf_s(fp_inicialization_data, "%lld %lld %lld %lld %lld %lld %lld %lld\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#else
			fprintf_s(fp_inicialization_data, "%d %d %d %d %d %d %d %d\n", inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8);
#endif
			}

		fclose(fp_inicialization_data);
	}
}


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
void exporttecplotxy360T_3D_part2_apparat_hot(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, integer ikey)
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

	FILE *fp;
	FILE *fp1; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal* temp_shadow = NULL;
	if (ionly_solid_visible == 1) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = NULL;
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
	switch (ikey) {
		case 0 : err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
		case 1 : err = fopen_s(&fp, "ALICEFLOW0_27_temp_apparat_hot.PLT", "wb");  break; // То что нужно для отчёта Алексею.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;  
	}
	if ((err) != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		char c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
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
					inode1 = database.nvtxcell[0][i] - 1;
					inode2 = database.nvtxcell[1][i] - 1;
					inode3 = database.nvtxcell[2][i] - 1;
					inode4 = database.nvtxcell[3][i] - 1;
					inode5 = database.nvtxcell[4][i] - 1;
					inode6 = database.nvtxcell[5][i] - 1;
					inode7 = database.nvtxcell[6][i] - 1;
					inode8 = database.nvtxcell[7][i] - 1;

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
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
						for (i = 0; i < database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = database.nvtxcell[0][i] - 1;
							inode2 = database.nvtxcell[1][i] - 1;
							inode3 = database.nvtxcell[2][i] - 1;
							inode4 = database.nvtxcell[3][i] - 1;
							inode5 = database.nvtxcell[4][i] - 1;
							inode6 = database.nvtxcell[5][i] - 1;
							inode7 = database.nvtxcell[6][i] - 1;
							inode8 = database.nvtxcell[7][i] - 1;
							
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
				if (bextendedprint) {
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				}
				else {
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
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
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.6f ", database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.6f ", database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.6f ", database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					else {
						// запись x
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.16f ", database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.16f ", database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.16f ", database.z[i]);
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

		}

			doublereal *Tx = NULL;
			doublereal *Ty = NULL;
			doublereal *Tz = NULL;
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
			if (Tx != NULL) {
				delete[] Tx;
			}
			if (Ty != NULL) {
				delete[] Ty;
			}
			if (Tz != NULL) {
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
			for (i = 0; i < database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = database.nvtxcell[0][i] - 1;
					inode2 = database.nvtxcell[1][i] - 1;
					inode3 = database.nvtxcell[2][i] - 1;
					inode4 = database.nvtxcell[3][i] - 1;
					inode5 = database.nvtxcell[4][i] - 1;
					inode6 = database.nvtxcell[5][i] - 1;
					inode7 = database.nvtxcell[6][i] - 1;
					inode8 = database.nvtxcell[7][i] - 1;

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
						//fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

						
					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

									}
						}
						else if ( ((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

								}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

									}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
					    {
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

								}
						}
					}
					else {
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
							// fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							
						}

					}
				}
			}
		}
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != NULL) {
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

	if (temp_shadow != NULL) {
		delete[] temp_shadow;
		temp_shadow = NULL;
	}
	//total_deformation
	
	if (total_deformation_shadow != NULL) {
		for (integer j_6 = 0; j_6 < 4; j_6++) {
			if (total_deformation_shadow[j_6] != NULL) {
				delete[] total_deformation_shadow[j_6];
				total_deformation_shadow[j_6] = NULL;
			}
		}
		delete[] total_deformation_shadow;
		total_deformation_shadow = NULL;
	}
	

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);



}

// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, integer ikey)
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

	FILE *fp;
	FILE *fp1; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal* temp_shadow = NULL;
	if (ionly_solid_visible == 1) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = NULL;
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
	switch (ikey) {
	case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp.PLT", "wb");  break; // То что нужно для отчёта Алексею.
	default: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	}
	if ((err) != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		char c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
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
					inode1 = database.nvtxcell[0][i] - 1;
					inode2 = database.nvtxcell[1][i] - 1;
					inode3 = database.nvtxcell[2][i] - 1;
					inode4 = database.nvtxcell[3][i] - 1;
					inode5 = database.nvtxcell[4][i] - 1;
					inode6 = database.nvtxcell[5][i] - 1;
					inode7 = database.nvtxcell[6][i] - 1;
					inode8 = database.nvtxcell[7][i] - 1;

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
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
						for (i = 0; i < database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = database.nvtxcell[0][i] - 1;
							inode2 = database.nvtxcell[1][i] - 1;
							inode3 = database.nvtxcell[2][i] - 1;
							inode4 = database.nvtxcell[3][i] - 1;
							inode5 = database.nvtxcell[4][i] - 1;
							inode6 = database.nvtxcell[5][i] - 1;
							inode7 = database.nvtxcell[6][i] - 1;
							inode8 = database.nvtxcell[7][i] - 1;

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
				// Полный набор искомых величин и теплопередача и гидродинамика:
				if (bextendedprint) {
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				}
				else {
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
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
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.6f ", database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.6f ", database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.6f ", database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					else {
						// запись x
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.16f ", database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись y
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.16f ", database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// запись z
						for (i = 0; i < database.maxelm; i++) {
							fprintf(fp, "%+.16f ", database.z[i]);
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

		}

		doublereal *Tx = NULL;
		doublereal *Ty = NULL;
		doublereal *Tz = NULL;
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
		if (Tx != NULL) {
			delete[] Tx;
		}
		if (Ty != NULL) {
			delete[] Ty;
		}
		if (Tz != NULL) {
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
			for (i = 0; i < database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = database.nvtxcell[0][i] - 1;
					inode2 = database.nvtxcell[1][i] - 1;
					inode3 = database.nvtxcell[2][i] - 1;
					inode4 = database.nvtxcell[3][i] - 1;
					inode5 = database.nvtxcell[4][i] - 1;
					inode6 = database.nvtxcell[5][i] - 1;
					inode7 = database.nvtxcell[6][i] - 1;
					inode8 = database.nvtxcell[7][i] - 1;

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
						//fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
					}
					else {
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
							// fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}
		}
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != NULL) {
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

	if (temp_shadow != NULL) {
		delete[] temp_shadow;
		temp_shadow = NULL;
	}
	//total_deformation

	if (total_deformation_shadow != NULL) {
		for (integer j_6 = 0; j_6 < 4; j_6++) {
			if (total_deformation_shadow[j_6] != NULL) {
				delete[] total_deformation_shadow[j_6];
				total_deformation_shadow[j_6] = NULL;
			}
		}
		delete[] total_deformation_shadow;
		total_deformation_shadow = NULL;
	}


	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);



} // apparat_hot



// 10 января 2016 . Заметка : надо сделать запись истинно бинарного файла, чтобы он быстрее открывался техплотом,
// а то при записи в текстовом режиме время открытия файла техплотом соизмеримо со временем вычисления. 
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2. Для анимации. Если inumbercadr==0 то первый кадр.
void exporttecplotxy360T_3D_part2_ianimation_series(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, integer ikey,
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

	FILE *fp;
	FILE *fp1; // часть 1 или 3
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	doublereal* temp_shadow = NULL;
	if (ionly_solid_visible == 1) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = NULL;
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
	
	if ((err) != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		char c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
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
					inode1 = database.nvtxcell[0][i] - 1;
					inode2 = database.nvtxcell[1][i] - 1;
					inode3 = database.nvtxcell[2][i] - 1;
					inode4 = database.nvtxcell[3][i] - 1;
					inode5 = database.nvtxcell[4][i] - 1;
					inode6 = database.nvtxcell[5][i] - 1;
					inode7 = database.nvtxcell[6][i] - 1;
					inode8 = database.nvtxcell[7][i] - 1;

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
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
						for (i = 0; i < database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = database.nvtxcell[0][i] - 1;
							inode2 = database.nvtxcell[1][i] - 1;
							inode3 = database.nvtxcell[2][i] - 1;
							inode4 = database.nvtxcell[3][i] - 1;
							inode5 = database.nvtxcell[4][i] - 1;
							inode6 = database.nvtxcell[5][i] - 1;
							inode7 = database.nvtxcell[6][i] - 1;
							inode8 = database.nvtxcell[7][i] - 1;

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
							for (i = 0; i < database.maxelm; i++) {
								fprintf(fp, "%+.6f ", database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись y
							for (i = 0; i < database.maxelm; i++) {
								fprintf(fp, "%+.6f ", database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись z
							for (i = 0; i < database.maxelm; i++) {
								fprintf(fp, "%+.6f ", database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						else {
							// запись x
							for (i = 0; i < database.maxelm; i++) {
								fprintf(fp, "%+.16f ", database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись y
							for (i = 0; i < database.maxelm; i++) {
								fprintf(fp, "%+.16f ", database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// запись z
							for (i = 0; i < database.maxelm; i++) {
								fprintf(fp, "%+.16f ", database.z[i]);
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

			}
	

		doublereal *Tx = NULL;
		doublereal *Ty = NULL;
		doublereal *Tz = NULL;
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
		if (Tx != NULL) {
			delete[] Tx;
		}
		if (Ty != NULL) {
			delete[] Ty;
		}
		if (Tz != NULL) {
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
			for (i = 0; i < database.ncell; i++) {
				if (bsolid_static_only) {
					//printf("Only solid ok\n");
					//getchar();

					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = database.nvtxcell[0][i] - 1;
					inode2 = database.nvtxcell[1][i] - 1;
					inode3 = database.nvtxcell[2][i] - 1;
					inode4 = database.nvtxcell[3][i] - 1;
					inode5 = database.nvtxcell[4][i] - 1;
					inode6 = database.nvtxcell[5][i] - 1;
					inode7 = database.nvtxcell[6][i] - 1;
					inode8 = database.nvtxcell[7][i] - 1;

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
						//fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
						fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и только как и раньше.
						//fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
								fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

							}
						}
					}
					else {
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;

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
							// fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
							fprintf(fp, "%lld %lld %lld %lld %lld %lld %lld %lld\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
							// Визуализация твёрдого тела и жидкости.
							// fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i], database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}
		}
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != NULL) {
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

	if (temp_shadow != NULL) {
		delete[] temp_shadow;
		temp_shadow = NULL;
	}
	//total_deformation

	if (total_deformation_shadow != NULL) {
		for (integer j_6 = 0; j_6 < 4; j_6++) {
			if (total_deformation_shadow[j_6] != NULL) {
				delete[] total_deformation_shadow[j_6];
				total_deformation_shadow[j_6] = NULL;
			}
		}
		delete[] total_deformation_shadow;
		total_deformation_shadow = NULL;
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
void exporttecplotxy360T_3D_part2amg(doublereal* u, bool bextendedprint, integer imove)
{
	integer ianimate = 0;
	integer flow_interior_count = 1;
	// imove 0 или 1 для нумерации с нуля или с единицы.

	// ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage = false;

	FILE *fp;
	FILE *fp1; // часть 1 или 3
	errno_t err;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
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
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		char c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
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
					//	fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", database.maxelm + t.maxbound, ncell);
#else
					//	fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", database.maxelm + t.maxbound, ncell);
#endif

				}
				else {

#if doubleintprecision == 1
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", database.maxelm, database.ncell);
#else
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", database.maxelm, database.ncell);
#endif

						}

				if (bvery_big_memory) {
					// extended printeger не предусмотрено.

					// запись x
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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
						for (i = 0; i < database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = database.nvtxcell[0][i] - 1;
							inode2 = database.nvtxcell[1][i] - 1;
							inode3 = database.nvtxcell[2][i] - 1;
							inode4 = database.nvtxcell[3][i] - 1;
							inode5 = database.nvtxcell[4][i] - 1;
							inode6 = database.nvtxcell[5][i] - 1;
							inode7 = database.nvtxcell[6][i] - 1;
							inode8 = database.nvtxcell[7][i] - 1;
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
					//fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", database.maxelm + t.maxbound, ncell_shadow);
#else
					//fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", database.maxelm + t.maxbound, ncell_shadow);
#endif

						}
				else {

#if doubleintprecision == 1
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", database.maxelm, database.ncell);
#else
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", database.maxelm, database.ncell);
#endif

					
				}
				if (bvery_big_memory) {
					// запись x
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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
			for (i = 0; i < database.maxelm; i++) {
				if (database.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
					fprintf(fp, "%+.16f ", u[database.ptr[0][i]+imove]); // my variable field
					avg_stat += u[database.ptr[0][i] + imove];
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			avg_stat /= database.maxelm;
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
		doublereal *Tx = NULL;
		doublereal *Ty = NULL;
		doublereal *Tz = NULL;
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
		if (Tx != NULL) {
			delete[] Tx;
		}
		if (Ty != NULL) {
			delete[] Ty;
		}
		if (Tz != NULL) {
			delete[] Tz;
		}
		*/


		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < database.ncell; i++) {
				if (bsolid_static_only) {

#if doubleintprecision == 1
					// Визуализация твёрдого тела и только как и раньше.
					fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
					// Визуализация твёрдого тела и только как и раньше.
					fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

						}
				else {
					if ((ionly_solid_visible == 1) && (flow_interior > 0))
					{
						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i] - 1;
						inode2 = database.nvtxcell[1][i] - 1;
						inode3 = database.nvtxcell[2][i] - 1;
						inode4 = database.nvtxcell[3][i] - 1;
						inode5 = database.nvtxcell[4][i] - 1;
						inode6 = database.nvtxcell[5][i] - 1;
						inode7 = database.nvtxcell[6][i] - 1;
						inode8 = database.nvtxcell[7][i] - 1;
						// визуализация только твёрдого тела.
						// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
						if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {

#if doubleintprecision == 1
							fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
							fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
							fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

								}
					}
					else {

#if doubleintprecision == 1
						// Визуализация твёрдого тела и жидкости.
						fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
						// Визуализация твёрдого тела и жидкости.
						fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

						}
				}
			}
		}
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != NULL) {
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

	FILE *fp;
	FILE *fp1; // часть 1 или 3
	errno_t err;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
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
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		char c; // читаемый символ
		integer ivarexport = 1; // по умолчанию только поле температур:
		integer i = 0; // счётчик цикла

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
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
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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
					for (i = 0; i < database.ncell; i++) {

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = database.nvtxcell[0][i]-1;
						inode2 = database.nvtxcell[1][i]-1;
						inode3 = database.nvtxcell[2][i]-1;
						inode4 = database.nvtxcell[3][i]-1;
						inode5 = database.nvtxcell[4][i]-1;
						inode6 = database.nvtxcell[5][i]-1;
						inode7 = database.nvtxcell[6][i]-1;
						inode8 = database.nvtxcell[7][i]-1;
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
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// запись z
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.z[i]);
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

		}

		doublereal *Tx = NULL;
		doublereal *Ty = NULL;
		doublereal *Tz = NULL;
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
		if (Tx != NULL) {
			delete[] Tx;
		}
		if (Ty != NULL) {
			delete[] Ty;
		}
		if (Tz != NULL) {
			delete[] Tz;
		}



		if (bvery_big_memory) {
			// запись информации о разностной сетке
			for (i = 0; i < database.ncell; i++) {
				if ((ionly_solid_visible == 1) && (flow_interior>0))
				{
					integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
					inode1 = database.nvtxcell[0][i]-1;
					inode2 = database.nvtxcell[1][i]-1;
					inode3 = database.nvtxcell[2][i]-1;
					inode4 = database.nvtxcell[3][i]-1;
					inode5 = database.nvtxcell[4][i]-1;
					inode6 = database.nvtxcell[5][i]-1;
					inode7 = database.nvtxcell[6][i]-1;
					inode8 = database.nvtxcell[7][i]-1;
					// визуализация только твёрдого тела.
					// идентификатор ==-1 говорит о том что узел принадлежит только твёрдому телу и не имеет жидкого представителя.
					if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
						(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {

#if doubleintprecision == 1
						fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
						fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

						}
				}
				else {


#if doubleintprecision == 1
					// Визуализация твёрдого тела и жидкости.
					fprintf(fp, "%lld %lld %lld %lld ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
					fprintf(fp, "%lld %lld %lld %lld\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#else
					// Визуализация твёрдого тела и жидкости.
					fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
#endif

					}
			}
		}
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				//getchar();
				system("pause");
				//exit(1);

			}
			else {

				if (fp1 != NULL) {
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
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, "xyplot.txt", "w");

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

	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, "xyplotT.txt", "w");

	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else {
		if (fp != NULL)
		{

			// внимание ! требуется указать точку через которую будет проходить линия и плоскость которой данная линия будет перпендикулярна.
			TOCHKA p;
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
			switch (iplane) {
			case XY: while (t.sosedi[BSIDE][iPf].iNODE1 < t.maxelm) iPf = t.sosedi[BSIDE][iPf].iNODE1; break;
			case XZ: while (t.sosedi[SSIDE][iPf].iNODE1 < t.maxelm) iPf = t.sosedi[SSIDE][iPf].iNODE1; break;
			case YZ: while (t.sosedi[WSIDE][iPf].iNODE1 < t.maxelm) iPf = t.sosedi[WSIDE][iPf].iNODE1; break;
			}

			integer G;
			switch (iplane) {
			case XY: G = TSIDE;  break;
			case XZ: G = NSIDE;  break;
			case YZ: G = ESIDE;  break;
			}


			fprintf(fp, "position,\ttemperature,\ttemperature_avg\n");
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
			volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
			center_cord3D(iPf, t.nvtx, t.pa, p,100); // вычисление координат центра КО.
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

void animationtecplot360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, char* title, bool btitle, integer iVar)
{
    // title - заголовок зоны.
	// btitle - печатать ли общий заголовок

	FILE *fp=NULL;
    FILE *fp1; // часть 1 или 3
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
	if (btitle) {
		// мы стираем предыдущие кадры и переходим к новой анимации:
		switch (iVar) {
	      case TEMP : err = fopen_s( &fp, "ALICEFLOW0_07_animation_temp.PLT",  "w");
		        break;
	      case SPEED : err = fopen_s( &fp, "ALICEFLOW0_07_animation_speed.PLT",  "w");
		         break;
	      case PRESS : err = fopen_s( &fp, "ALICEFLOW0_07_animation_press.PLT",  "w");
		         break;
	      case PAM :  err = fopen_s( &fp, "ALICEFLOW0_07_animation_pam.PLT",  "w");
		        break;
	   }
	}
	else {
		// мы добавляем следующие анимационные кадры.
	    switch (iVar) {
	      case TEMP : err = fopen_s( &fp, "ALICEFLOW0_07_animation_temp.PLT",  "a");
		        break;
	      case SPEED : err = fopen_s( &fp, "ALICEFLOW0_07_animation_speed.PLT",  "a");
		         break;
	      case PRESS : err = fopen_s( &fp, "ALICEFLOW0_07_animation_press.PLT",  "a");
		         break;
	      case PAM :  err = fopen_s( &fp, "ALICEFLOW0_07_animation_pam.PLT",  "a");
		        break;
	   }
	}
	
	if (err != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {
        
        char c; // читаемый символ
		integer ivarexport=1; // по умолчанию только поле температур:
		integer i=0; // счётчик цикла

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
		    printf("Open File temp part1 Error\n");
		    //getchar();
			system("pause");

	    }
	    else {
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
				
		

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

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
                

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

						
			}
            fclose(fp1); // закрытие файла
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

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
		    printf("Open File temp part3 Error\n");
		    //getchar();
			system("pause");

	    }
	    else {
			if (fp1 != NULL) {
				// копирование третьей части в итоговый файл
				while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				fclose(fp1); // закрытие файла
				//printf("export tecplot part1 is successfully reading and written...OK.\n");
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
void animationtecplot360T_3D_part2all(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, char* title, bool btitle)
{
	// величины которые следует анимировать.
	bool bTEMP=true;
	bool bSpeed=true;
	bool bPressure=true;
	bool bPAM=true;

	// Непосредственно запись анимации в разные файлы.
	if (bTEMP) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, TEMP);
	}

	if (bSpeed) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, SPEED);
	}

	if (bPressure) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, PRESS);
	}

	if (bPAM) {
        animationtecplot360T_3D_part2(maxelm, ncell, f, t, flow_interior_count, title, btitle, PAM);
	}
	
}

#endif

// Файл my_export_tecplot3.c передача результатов
// моделирования в программу tecplot360

#ifndef MY_EXPORT_TECPLOT3_C
#define MY_EXPORT_TECPLOT3_C 1

//#include "windows.h" // для функции WinExec
#include <string.h>

// проверка построенной сетки
// экспорт результата расчёта в программу tecplot360
void exporttecplotxy360_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal** potent, doublereal **rhie_chow)
{
	FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "ALICEFLOW0_03.PLT", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

		// запись имён переменных
		fprintf(fp, "VARIABLES = x, y, z, Vx, Vy, Vz, Mag, Pressure , Normal , PAm, Zero\n");

		// запись информации о зонах
		//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
        //if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
        fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);

		integer i=0; // счётчики 
		integer j=0; // цикла for

		// запись x
	    for (i=0; i<maxelm; i++) {	
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x));
			if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// запись y
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].y+pa[nvtx[2][i]-1].y));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

        // запись z
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		
		// запись горизонтальной Vx скорости
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[VX][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

        // запись вертикальной Vy скорости
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[VY][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

        // запись Vz скорости
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[VZ][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");


		// запись модуля скорости
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", sqrt(potent[VX][i]*potent[VX][i]+potent[VY][i]*potent[VY][i]+potent[VZ][i]*potent[VZ][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// запись давления
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[PRESS][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");
        
		// запись Normal
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", rhie_chow[0][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// запись PAm
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ",potent[PAM][i] ); // rhie_chow[1][i]
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// запись Rhi-Chow
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", rhie_chow[2][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// запись информации о разностной сетке
		for (i=0;i<ncell; i++) {
			fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
            fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
		}

		fclose(fp); // закрытие файла
		printf("file is successfully written...OK.\n");
       // WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
		  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);
	}
} // exporttecplotxy360_3D

// проверка построенной сетки
// экспорт результата расчёта в программу tecplot360
void exporttecplotxy360T_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal* potent)
{
	FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "ALICEFLOW0_03.PLT", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

		// запись имён переменных
		fprintf(fp, "VARIABLES = x, y, z, Temp\n");

		// запись информации о зонах
        fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);

		integer i=0; // счётчики 
		integer j=0; // цикла for

		// запись x
	    for (i=0; i<maxelm; i++) {	
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x));
			if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// запись y
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].y+pa[nvtx[2][i]-1].y));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

        // запись z
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		
		// запись температуры
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// запись информации о разностной сетке
		for (i=0;i<ncell; i++) {
			fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
            fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
		}

		fclose(fp); // закрытие файла
		printf("file is successfully written...OK.\n");
       // WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
		  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);
	}
} // exporttecplotxy360T_3D

// Трёхэтапная запись файла позволяет сэкономить 19N оперативной памяти.
// Далее везде применяется трёхэтампное формирование выходного файла.

// проверка построенной сетки
// экспорт результата расчёта в программу tecplot360
// части 1 и 3.
void exporttecplotxy360T_3D_part1and3(integer maxelm, integer maxbound, bool bextendedprint, integer ncell, 
									  integer** nvtx, integer** nvtxcell, TOCHKA* pa,
									  BOUND* border_neighbor, integer ivarexport)
{
	// расширенная печать
	// При расширенной печати мы печатаем также и граничные узлы.
	// bextendedprint=true; расширенная печать. 

    // ivarexport == 1 печатается только поле температур,
	// ivarexport == 2 печатается только гидродинамика,
	// ivarexport == 3 печатается и поле температур и гидродинамика.

	FILE *fp;
	errno_t err;
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	// запись частей 1 и 3
	if ((err = fopen_s( &fp, "ALICEFLOW0_06_temp_part1.txt", "w")) != 0) {
		printf("Create File temp part1 Error\n");
		getchar();

	}
	else {
		// запись заголовка
		/*fprintf(fp, "TITLE = \"ALICEFLOW0_06\"\n");

		// запись имён переменных
		switch (ivarexport) {
			case 1: fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n"); break; // печатается только поле температур
			case 2: fprintf(fp, "VARIABLES = x, y, z, Speed, Pressure, Vx, Vy, Vz, Rho, Mu, Mut\n"); break; // печатается только гидродинамика
			case 3: fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, Speed, Pressure, Vx, Vy, Vz, Rho, Mu, Mut\n"); break; // печатается и температура и гидродинамика
			default: printf("Error export tecplot. Nonselected exporting variables...\n"); getchar();
		}

		// запись информации о зонах
        fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
		*/

		integer i=0; // счётчики 
		integer j=0; // цикла for
		doublereal m2cm=100.0;

		// запись x
	    for (i=0; i<maxelm; i++) {	
			fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x)*m2cm);
			if (i%10==0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i=0; i<maxbound; i++) {
				switch (border_neighbor[i].Norm) {// определим внутреннюю нормаль к границе
				case E: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI]-1].x*m2cm);
					     break;
				case W: fprintf(fp, "%+.16f ", pa[nvtx[1][border_neighbor[i].iI]-1].x*m2cm);
					     break;
				case N: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].x+pa[nvtx[1][border_neighbor[i].iI]-1].x)*m2cm);
					     break;
				case S: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].x+pa[nvtx[1][border_neighbor[i].iI]-1].x)*m2cm);
					     break;
				case T: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].x+pa[nvtx[1][border_neighbor[i].iI]-1].x)*m2cm);
					     break;
				case B: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].x+pa[nvtx[1][border_neighbor[i].iI]-1].x)*m2cm);
					     break;
				}
			    if (i%10==0) fprintf(fp, "\n");
			}
		}
			
		fprintf(fp, "\n");
          
		// запись y
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i]-1].y+pa[nvtx[2][i]-1].y)*m2cm);
            if (i%10==0) fprintf(fp, "\n");
		}
		
		if (bextendedprint) {
			for (i=0; i<maxbound; i++) {
				switch (border_neighbor[i].Norm) {// определим внутреннюю нормаль к границе
				case E: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].y+pa[nvtx[2][border_neighbor[i].iI]-1].y)*m2cm);
					     break;
				case W: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].y+pa[nvtx[2][border_neighbor[i].iI]-1].y)*m2cm);
					     break;
				case N: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI]-1].y*m2cm);
					     break;
				case S: fprintf(fp, "%+.16f ", pa[nvtx[2][border_neighbor[i].iI]-1].y*m2cm);
					     break;
				case T: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].y+pa[nvtx[2][border_neighbor[i].iI]-1].y)*m2cm);
					     break;
				case B: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].y+pa[nvtx[2][border_neighbor[i].iI]-1].y)*m2cm);
					     break;
				}
			    if (i%10==0) fprintf(fp, "\n");
			}
		}

        fprintf(fp, "\n");

        // запись z
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z)*m2cm);
            if (i%10==0) fprintf(fp, "\n");
		}
		
		if (bextendedprint) {
			for (i=0; i<maxbound; i++) {
				switch (border_neighbor[i].Norm) {// определим внутреннюю нормаль к границе
				case E: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].z+pa[nvtx[4][border_neighbor[i].iI]-1].z)*m2cm);
					     break;
				case W: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].z+pa[nvtx[4][border_neighbor[i].iI]-1].z)*m2cm);
					     break;
				case N: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].z+pa[nvtx[4][border_neighbor[i].iI]-1].z)*m2cm);
					     break;
				case S: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI]-1].z+pa[nvtx[4][border_neighbor[i].iI]-1].z)*m2cm);
					     break;
				case T: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI]-1].z*m2cm);
					     break;
				case B: fprintf(fp, "%+.16f ", pa[nvtx[4][border_neighbor[i].iI]-1].z*m2cm);
					     break;
				}
			    if (i%10==0) fprintf(fp, "\n");
			}
		}

        fprintf(fp, "\n");

		fclose(fp); // закрытие файла
		printf("export tecplot temperature part1 is successfully written...OK.\n");
       
	}

    if ((err = fopen_s( &fp, "ALICEFLOW0_06_temp_part3.txt", "w")) != 0) {
		printf("Create File temp part3 Error\n");
		getchar();
	}
	else {

		integer i=0;
		// запись информации о разностной сетке
		for (i=0;i<ncell; i++) {
			fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
            fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
		}
        
		fclose(fp); // закрытие файла
		printf("export tecplot temperature part3 is successfully written...OK.\n");
	}

} // exporttecplotxy360T_3D_part1and3


// построение графика температуры вдоль линии 
void xyplot_temp(TEMPER &t, doublereal* tempfiltr) {

	// tempfiltr - передаваемая однократно фильтрованная температура.
	// имя создаваемого файла xyplotT.txt.

	FILE *fp;
	errno_t err;
	if ((err = fopen_s( &fp, "xyplotT.txt",  "w")) != 0) {
		printf("Create File xyplot Error\n");
		getchar();

	}
	else {

		// внимание ! требуется указать точку через которую будет проходить линия и плоскость которой данная линия будет перпендикулярна.
		TOCHKA p;
		doublereal epsilon=1e30;
		doublereal dist;
		doublereal x=0.0, y=0.0, z=10.0e-6; // точка через которую проходит линия
		integer iPf=0;
		integer iplane=XY; // плоскость перпендикулярная линии.
		
		for (integer iP=0; iP<t.maxelm; iP++) {
			center_cord3D(iP, t.nvtx, t.pa, p); // вычисление координат центра КО.
			dist=sqrt(fabs(x-p.x)*fabs(x-p.x)+fabs(y-p.y)*fabs(y-p.y)+fabs(z-p.z)*fabs(z-p.z));
			if (dist<epsilon) {
				epsilon=dist;
				iPf=iP;
			}
		}
		
		// перемотка в начало 
		switch (iplane) {
		  case XY: while (t.neighbors_for_the_internal_node[B][iPf].iNODE1<t.maxelm) iPf=t.neighbors_for_the_internal_node[B][iPf].iNODE1; break;
		  case XZ: while (t.neighbors_for_the_internal_node[S][iPf].iNODE1<t.maxelm) iPf=t.neighbors_for_the_internal_node[S][iPf].iNODE1; break;
		  case YZ: while (t.neighbors_for_the_internal_node[W][iPf].iNODE1<t.maxelm) iPf=t.neighbors_for_the_internal_node[W][iPf].iNODE1; break;
		}

		integer G;
		switch (iplane) {
		  case XY: G=T;  break;
		  case XZ: G=N;  break;
		  case YZ: G=E;  break;
		}
		
		
		fprintf(fp, "position,\ttemperature,\ttemperature_avg\n");
		doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	    volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
        center_cord3D(iPf, t.nvtx, t.pa, p); // вычисление координат центра КО.
		switch (iplane) {
		  case XY: fprintf(fp, "%+.16f %+.16f %+.16f\n",
			  p.z - 0.5*dz, t.potent[t.neighbors_for_the_internal_node[B][iPf].iNODE1],
			  tempfiltr[t.neighbors_for_the_internal_node[B][iPf].iNODE1]);
			              break;
		  case XZ:  fprintf(fp,"%+.16f %+.16f %+.16f\n",
			  p.y - 0.5*dy, t.potent[t.neighbors_for_the_internal_node[S][iPf].iNODE1],
			  tempfiltr[t.neighbors_for_the_internal_node[S][iPf].iNODE1]);
			              break;
		  case YZ:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
			  p.x - 0.5*dx, t.potent[t.neighbors_for_the_internal_node[W][iPf].iNODE1],
			  tempfiltr[t.neighbors_for_the_internal_node[W][iPf].iNODE1]);
			              break;
		}
		switch (iplane) {
		  case XY: while (iPf<t.maxelm) {
			        center_cord3D(iPf, t.nvtx,t.pa, p); 
			        fprintf(fp, "%+.16f %+.16f %+.16f\n", 
						  p.z, t.potent[iPf],
					      tempfiltr[iPf]);
					if (t.neighbors_for_the_internal_node[T][iPf].iNODE1 >= t.maxelm) {
						  volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f\n",
							  p.z + 0.5*dz, t.potent[t.neighbors_for_the_internal_node[T][iPf].iNODE1],
							  tempfiltr[t.neighbors_for_the_internal_node[T][iPf].iNODE1]);
					}
					iPf = t.neighbors_for_the_internal_node[T][iPf].iNODE1;
					} break;
		  case XZ: while (iPf<t.maxelm) {
			        center_cord3D(iPf, t.nvtx, t.pa, p); 
			        fprintf(fp, "%+.16f %+.16f %+.16f\n", 
						  p.y, t.potent[iPf],
					      tempfiltr[iPf]);
					if (t.neighbors_for_the_internal_node[N][iPf].iNODE1 >= t.maxelm) {
						  volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f\n",
							  p.y + 0.5*dy, t.potent[t.neighbors_for_the_internal_node[N][iPf].iNODE1],
							  tempfiltr[t.neighbors_for_the_internal_node[N][iPf].iNODE1]);
					}
					/*
					// Узнаём последовательность узлов для отладки.
					printf("iPf=%d\n",iPf);
					if (fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1>=fglobal[ifi].maxelm) {
						printf("iPffinish=%d\n",fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1);
						getchar();
					}
					*/
					iPf = t.neighbors_for_the_internal_node[N][iPf].iNODE1;
					} break;
		  case YZ: while (iPf<t.maxelm) {
			        center_cord3D(iPf, t.nvtx, t.pa, p); 
			        fprintf(fp, "%+.16f %+.16f %+.16f\n", 
						  p.y, t.potent[iPf],
					      tempfiltr[iPf]);
					if (t.neighbors_for_the_internal_node[E][iPf].iNODE1 >= t.maxelm) {
						  volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f\n",
							  p.x + 0.5*dx, t.potent[t.neighbors_for_the_internal_node[E][iPf].iNODE1],
							  tempfiltr[t.neighbors_for_the_internal_node[E][iPf].iNODE1]);
					}
					iPf = t.neighbors_for_the_internal_node[E][iPf].iNODE1;
					} break;
		}
		fclose(fp);
	}
} // xyplot_temp


// проверка построенной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint)
{
    // ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage=false;

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

	
	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	if ((err = fopen_s( &fp, "ALICEFLOW0_07_temp.PLT",  "w")) != 0) {
		printf("Create File temp Error\n");
		getchar();

	}
	else {
        
        char c; // читаемый символ
		integer ivarexport=1; // по умолчанию только поле температур:
		integer i=0; // счётчик цикла

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
		    printf("Open File temp part1 Error\n");
		    getchar();

	    }
	    else {
			// копирование первой части в итоговый файл
			// Особенность: иногда необходимо изменить вторую строку в файле:
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
		        /// fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, heat_flux_x, heat_flux_y, heat_flux_z\n");  // печатается только поле температур
				  fprintf(fp, "VARIABLES = x_cm, y_cm, z_cm, electric_potencial_V, EPS, Ex_V/cm, Ey_V/cm, Ez_V/cm, Emag_V/cm\n");
			
				 // запись информации о зонах
				 if (bextendedprint) {
					 fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm+t.maxbound, ncell);
				 }
				 else {
                    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				 }
		

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);
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
					
				// запись информации о зонах
				if (bextendedprint) {
                    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm+t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);
			}
            fclose(fp1); // закрытие файла
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// запись второй части
        
        // Запись поля температур производится всегда.


		// запись температуры
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.potent[i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i=maxelm;i<maxelm+t.maxbound; i++) {
			    fprintf(fp, "%+.16f ", t.potent[i]);
                if ((i+maxelm)%10==0) fprintf(fp, "\n");
			}
		}
		

		fprintf(fp, "\n");

		// Lam
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%+.16f ", t.prop[LAM][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		if (bextendedprint) {
			for (i=0;i<t.maxbound; i++) {
			     fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
                 if ((i+maxelm)%10==0) fprintf(fp, "\n");
		    }
		}

		fprintf(fp, "\n");

		if (ivarexport==1) {

		/*
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
					t.neighbors_for_the_internal_node, t.maxelm, false, 
					t.border_neighbor,  Tx, Ty, Tz);
			}

			for (i=0; i<t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, true, 
					t.border_neighbor,Tx, Ty, Tz);
			}
			

			// Сохранение в файл.

			// Heat Flux X
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Tx[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Tx[i+maxelm]);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Heat Flux Y
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Ty[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Ty[i+maxelm]);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Heat Flux Z
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Tz[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Tz[i+maxelm]);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Освобождение оперативной памяти.
			delete Tx;
			delete Ty;
			delete Tz;

			*/

			doublereal m2cm=100.0;

		doublereal *Tx=NULL;
            doublereal *Ty=NULL;
            doublereal *Tz=NULL;
			doublereal *tempfiltr=NULL;
			Tx=new doublereal[t.maxelm+t.maxbound];
			Ty=new doublereal[t.maxelm+t.maxbound];
			Tz=new doublereal[t.maxelm+t.maxbound];
			tempfiltr=new doublereal[t.maxelm+t.maxbound];

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
					t.neighbors_for_the_internal_node, t.maxelm, false, 
					t.border_neighbor,  Tx, Ty, Tz);
			}

			for (i=0; i<t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, true, 
					t.border_neighbor,Tx, Ty, Tz);
			}
			

			// Сохранение в файл.

			// EX
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -Tx[i]/m2cm);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -Tx[i+maxelm]/m2cm);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// EY
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -Ty[i]/m2cm);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -Ty[i+maxelm]/m2cm);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// EZ
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -Tz[i]/m2cm);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -Tz[i+maxelm]/m2cm);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Emag
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", sqrt(Tx[i]*Tx[i]+Ty[i]*Ty[i]+Tz[i]*Tz[i])/m2cm);
                if (i%10==0) fprintf(fp, "\n");
				tempfiltr[i]=sqrt(Tx[i]*Tx[i]+Ty[i]*Ty[i]+Tz[i]*Tz[i])/m2cm;
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", sqrt(Tx[i+maxelm]*Tx[i+maxelm]+Ty[i+maxelm]*Ty[i+maxelm]+Tz[i+maxelm]*Tz[i+maxelm])/m2cm);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
				   tempfiltr[i+maxelm]=sqrt(Tx[i+maxelm]*Tx[i+maxelm]+Ty[i+maxelm]*Ty[i+maxelm]+Tz[i+maxelm]*Tz[i+maxelm])/m2cm;
		       }
		    }

		    fprintf(fp, "\n");

			xyplot_temp(t, tempfiltr);

			// Освобождение оперативной памяти.
			delete Tx;
			delete Ty;
			delete Tz;
			delete tempfiltr;

			}


		// Запись гидродинамических величин если необходимо:
		if (ivarexport==3) {
			// Speed
            for (i=0;i<maxelm; i++) {
				if (t.ptr[1][i]>-1) {
					doublereal svx=f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VX][t.ptr[0][i]];
					doublereal svy=f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VY][t.ptr[0][i]];
					doublereal svz=f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]];
					fprintf(fp, "%+.16f ",sqrt(svx+svy+svz)); 
				} else fprintf(fp, "%+.16f ", 0.0);
                 if (i%10==0) fprintf(fp, "\n");
		    }

			
			if (bextendedprint) {
				// Speed
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
				    	doublereal svx=f[idfluid].potent[VX][i+maxelm]*f[idfluid].potent[VX][i+maxelm];
					    doublereal svy=f[idfluid].potent[VY][i+maxelm]*f[idfluid].potent[VY][i+maxelm];
					    doublereal svz=f[idfluid].potent[VZ][i+maxelm]*f[idfluid].potent[VZ][i+maxelm];
					    fprintf(fp, "%+.16f ",sqrt(svx+svy+svz)); 
                      if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}
			

		    fprintf(fp, "\n");

			// Pressure
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// Pressure
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i+maxelm]); // PRESSURE
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// PAM
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// PAM
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i+maxelm]); // PRESSURE
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// VX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VX][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// VX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[VX][i+maxelm]); // VX
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// VY
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VY][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// VY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[VY][i+maxelm]); // VY
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// VZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VZ][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// VZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[VZ][i+maxelm]); // VZ
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			
            // Rho
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VX][i]);
				} else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// Rho
				for (i=0;i<f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ",f[0].prop_b[RHO][i]);
					if ((i+maxelm)%10==0) fprintf(fp, "\n");
				}
			}

		    fprintf(fp, "\n");

			// Mu
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VX][i].ap);
				} else fprintf(fp, "%+.16f ", 0.0); 
				// Вязкость в твёрдом теле при визуализации положена равной нулю, хотя должна быть бесконечно большой.
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// Mu
				for (i=0;i<f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ",f[0].prop_b[MU][i]);
					if ((i+maxelm)%10==0) fprintf(fp, "\n");
				}
			}

		    fprintf(fp, "\n");

			// Mut // Турбулентная динамическая вязкость
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// MUT
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i+maxelm]); // MUT
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// для отладки график распределения нумерации контрольных объёмов.
			// или распределение расстояния до стенки Distance_Wall.
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   //fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime==ZEROEQMOD) || (f[t.ptr[1][i]].iflowregime==SMAGORINSKY)) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].rdistWall[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// Distance_Wall.
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
					if ((f[0].iflowregime==ZEROEQMOD) || (f[0].iflowregime==SMAGORINSKY)) {
                       fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i+maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			

			// Curl // Завихрённость - модуль ротора скорости.
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// Curl
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i+maxelm]); // Curl
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// Частные производные от компонент скорости !!!.

			// GRADXVX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADXVX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i+maxelm]); // GRADXVX
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADYVX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADYVX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i+maxelm]); // GRADYVX
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADZVX
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADZVX
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i+maxelm]); // GRADZVX
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADXVY
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADXVY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i+maxelm]); // GRADXVY
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADYVY
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADYVY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i+maxelm]); // GRADYVY
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADZVY
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADZVY
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i+maxelm]); // GRADZVY
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADXVZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADXVZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i+maxelm]); // GRADXVZ
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADYVZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADYVZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i+maxelm]); // GRADYVZ
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

			// GRADZVZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]); 
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

			if (bextendedprint) {
				// GRADZVZ
				integer idfluid=0;
                for (i=0;i<f[idfluid].maxbound; i++) {
                       fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVZ][i+maxelm]); // GRADZVZ
                    if ((i+maxelm)%10==0) fprintf(fp, "\n");
		        }
			}

		    fprintf(fp, "\n");

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
					t.neighbors_for_the_internal_node, t.maxelm, false, 
					t.border_neighbor,  Tx, Ty, Tz);
			}

			for (i=0; i<t.maxelm; i++) {
				// Только граничные узлы.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, true, 
					t.border_neighbor,Tx, Ty, Tz);
			}
			

			// Сохранение в файл.

			// Heat Flux X
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Tx[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Tx[i+maxelm]);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Heat Flux Y
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Ty[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Ty[i+maxelm]);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Heat Flux Z
		    for (i=0;i<maxelm; i++) {
			    fprintf(fp, "%+.16f ", -t.prop[LAM][i]*Tz[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    if (bextendedprint) {
			   for (i=0;i<t.maxbound; i++) {
			       fprintf(fp, "%+.16f ", -t.prop_b[LAM][i]*Tz[i+maxelm]);
                   if ((i+maxelm)%10==0) fprintf(fp, "\n");
		       }
		    }

		    fprintf(fp, "\n");

			// Освобождение оперативной памяти.
			delete Tx;
			delete Ty;
			delete Tz;

		}

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
		    printf("Open File temp part3 Error\n");
		    getchar();

	    }
	    else {
			// копирование третьей части в итоговый файл
			while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);
            fclose(fp1); // закрытие файла
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		fclose(fp); // закрытие файла
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}

void xyplot( FLOW* &fglobal, integer flow_interior, TEMPER &t) {
	FILE *fp;
	errno_t err;
	if ((err = fopen_s( &fp, "xyplot.txt",  "w")) != 0) {
		printf("Create File xyplot Error\n");
		getchar();

	}
	else {

		TOCHKA p;
		doublereal epsilon=1e30;
		doublereal dist;
		doublereal x=0.0e-3, y=0.0e-3, z=0.0e-3; // точка через которую проходит линия
		integer ifi=0, iPf=0;
		integer iplane=XZ; // плоскость перпендикулярная линии.
		for (integer i=0; i<flow_interior; i++) {
			for (integer iP=0; iP<fglobal[i].maxelm; iP++) {
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p); // вычисление координат центра КО.
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
		case XY: while (fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1<fglobal[ifi].maxelm) iPf = fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1; break;
		case XZ: while (fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1<fglobal[ifi].maxelm) iPf = fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1; break;
		case YZ: while (fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1<fglobal[ifi].maxelm) iPf = fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1; break;
		}

		integer G;
		switch (iplane) {
		  case XY: G=T;  break;
		  case XZ: G=N;  break;
		  case YZ: G=E;  break;
		}
		
		
		fprintf(fp, "position,\tVx,\tVy,\tVz,\tPam,\tPress,\tFbuf,\tGRADPRESS,\tGRADPAM,\tmassfluxingran\n");
		doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	    volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p); // вычисление координат центра КО.
		switch (iplane) {
		  case XY: fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			  p.z - 0.5*dz, fglobal[ifi].potent[VX][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[VY][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[VZ][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[GRADZPRESS][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
			  fglobal[ifi].potent[GRADZPAM][fglobal[ifi].neighbors_for_the_internal_node[B][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
			              break;
		  case XZ:  fprintf(fp,"%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			  p.y - 0.5*dy, fglobal[ifi].potent[VX][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[VY][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[VZ][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[GRADYPRESS][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
			  fglobal[ifi].potent[GRADYPAM][fglobal[ifi].neighbors_for_the_internal_node[S][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
			              break;
		  case YZ:  fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			  p.x - 0.5*dx, fglobal[ifi].potent[VX][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[VY][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[VZ][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[GRADXPRESS][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
			  fglobal[ifi].potent[GRADXPAM][fglobal[ifi].neighbors_for_the_internal_node[W][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
			              break;
		}
		switch (iplane) {
		  case XY: while (iPf<fglobal[ifi].maxelm) {
			        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p); 
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
					if (fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1 >= fglobal[ifi].maxelm) {
						  volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
							  p.z + 0.5*dz, fglobal[ifi].potent[VX][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[VY][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[VZ][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[GRADZPRESS][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
							  fglobal[ifi].potent[GRADZPAM][fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
					}
					iPf = fglobal[ifi].neighbors_for_the_internal_node[T][iPf].iNODE1;
					} break;
		  case XZ: while (iPf<fglobal[ifi].maxelm) {
			        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p); 
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
					if (fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1 >= fglobal[ifi].maxelm) {
						  volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
							  p.y + 0.5*dy, fglobal[ifi].potent[VX][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[VY][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[VZ][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[GRADYPRESS][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
							  fglobal[ifi].potent[GRADYPAM][fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
					}
					/*
					// Узнаём последовательность узлов для отладки.
					printf("iPf=%d\n",iPf);
					if (fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1>=fglobal[ifi].maxelm) {
						printf("iPffinish=%d\n",fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1);
						getchar();
					}
					*/
					iPf = fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1;
					} break;
		  case YZ: while (iPf<fglobal[ifi].maxelm) {
			        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p); 
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
					if (fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1 >= fglobal[ifi].maxelm) {
						  volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
                          fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
							  p.x + 0.5*dx, fglobal[ifi].potent[VX][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[VY][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[VZ][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[GRADXPRESS][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
							  fglobal[ifi].potent[GRADXPAM][fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1],
						  fglobal[ifi].mf[iPf][G]);
					}
					iPf = fglobal[ifi].neighbors_for_the_internal_node[E][iPf].iNODE1;
					} break;
		}
		fclose(fp);
	}
}




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
	      case TEMP: err = fopen_s( &fp, "ALICEFLOW0_07_animation_temp.PLT",  "w");
		        break;
	      case SPEED: err = fopen_s( &fp, "ALICEFLOW0_07_animation_speed.PLT",  "w");
		         break;
	      case PRESS: err = fopen_s( &fp, "ALICEFLOW0_07_animation_press.PLT",  "w");
		         break;
	      case PAM:  err = fopen_s( &fp, "ALICEFLOW0_07_animation_pam.PLT",  "w");
		        break;
	   }
	}
	else {
		// мы добавляем следующие анимационные кадры.
	    switch (iVar) {
	      case TEMP: err = fopen_s( &fp, "ALICEFLOW0_07_animation_temp.PLT",  "a");
		        break;
	      case SPEED: err = fopen_s( &fp, "ALICEFLOW0_07_animation_speed.PLT",  "a");
		         break;
	      case PRESS: err = fopen_s( &fp, "ALICEFLOW0_07_animation_press.PLT",  "a");
		         break;
	      case PAM:  err = fopen_s( &fp, "ALICEFLOW0_07_animation_pam.PLT",  "a");
		        break;
	   }
	}
	
	if (err != 0) {
		printf("Create File temp Error\n");
		getchar();

	}
	else {
        
        char c; // читаемый символ
		integer ivarexport=1; // по умолчанию только поле температур:
		integer i=0; // счётчик цикла

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
		    printf("Open File temp part1 Error\n");
		    getchar();

	    }
	    else {
			// копирование первой части в итоговый файл
			// Особенность: иногда необходимо изменить вторую строку в файле:
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
	               case TEMP: fprintf(fp, "VARIABLES = x, y, z, Temp\n");  // печатается только поле температур
		                       break;
	               case SPEED:  fprintf(fp, "VARIABLES = x, y, z, Speed\n"); // печатается модуль скорости
		                       break;
	               case PRESS:  fprintf(fp, "VARIABLES = x, y, z, Press\n"); // печатается давление
		                       break;
	               case PAM:   fprintf(fp, "VARIABLES = x, y, z, PAM\n"); // печатается поправка давления
		                       break;
	            }
		         
			
				 // запись информации о зонах
                fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
		

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

			}
			else if (ivarexport==3) {
				
				if (btitle) {

				    // запись заголовка
		            fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

     				// запись имён переменных
				    switch (iVar) {
	                   case TEMP: fprintf(fp, "\nVARIABLES = x, y, z, Temp\n");  // печатается только поле температур
		                           break;
	                   case SPEED:  fprintf(fp, "\nVARIABLES = x, y, z, Speed\n"); // печатается модуль скорости
		                           break;
	                   case PRESS:  fprintf(fp, "\nVARIABLES = x, y, z, Press\n"); // печатается давление
		                           break;
	                   case PAM:   fprintf(fp, "\nVARIABLES = x, y, z, PAM\n"); // печатается поправка давления
		                           break;
	                }
				}
					
				// запись информации о зонах
                fprintf(fp, "ZONE T=\"%s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", title,  maxelm, ncell);

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
		    getchar();

	    }
	    else {
			// копирование третьей части в итоговый файл
			while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);
            fclose(fp1); // закрытие файла
            //printf("export tecplot part1 is successfully reading and written...OK.\n");
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

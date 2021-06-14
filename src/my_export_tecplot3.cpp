// ���� my_export_tecplot3.c �������� �����������
// ������������� � ��������� tecplot360

#pragma once
#ifndef MY_EXPORT_TECPLOT3_CPP
#define MY_EXPORT_TECPLOT3_CPP 1

//#include "windows.h" // ��� ������� WinExec
#include <string.h>
#include <stdlib.h> // for ������������ qsort
#include <algorithm> // sort (IntroSort)


// ��������� ������ ��� ������������������ ��������.
typedef struct TSORT_XY_PLOT_FLOW {

	doublereal x_argument;
	doublereal Vx;
	doublereal Vy;
	doublereal Vz;
	doublereal Speed;
	doublereal Pam;
	doublereal Press;
	doublereal Fbuf;
	doublereal GradPress; // ����� ����� ���������� �������.
	doublereal GradPam;  // ����� ����� ���������� �������.
	doublereal massflux_in_gran;  // ����� ����� ���������� �������.
								  // additional data:
								  // ����������� �� ��������� ��������:
	doublereal GradXVx;
	doublereal GradYVx;
	doublereal GradZVx;
	doublereal GradXVy;
	doublereal GradYVy;
	doublereal GradZVy;
	doublereal GradXVz;
	doublereal GradYVz;
	doublereal GradZVz;
	// ������������ � ��������� ��������.
	doublereal Curl;
	doublereal GradXPress;
	doublereal GradYPress;
	doublereal GradZPress;
	// ��������� ������������� �������.
	// ������ ��������-���������.
	doublereal Nusha;
	doublereal GradXNusha;
	doublereal GradYNusha;
	doublereal GradZNusha;
	// ������ SST K-Omega �������.
	doublereal Ke;
	doublereal Omega;
	doublereal GradXKe;
	doublereal GradYKe;
	doublereal GradZKe;
	doublereal GradXOmega;
	doublereal GradYOmega;
	doublereal GradZOmega;

	doublereal F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16;
	doublereal F17, F18, F19, F20, F21, F22, F23, F24, F25, F26, F27, F28, F29, F30, F31, F32;


} SORT_XY_PLOT_FLOW;

bool compare_XY_PLOT_FLOW(SORT_XY_PLOT_FLOW& Amat1, SORT_XY_PLOT_FLOW& Amat2) {
	return (Amat1.x_argument < Amat2.x_argument);
} // compare_XY_PLOT


// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
void exporttecplotxy360_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal** potent, doublereal **rhie_chow)
{
	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("ALICEFLOW0_03.PLT", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_03.PLT", "w");
#endif

	
	// �������� ����� ��� ������.
	if ((err) != 0) {
		printf("Create File Error\n");
		exit(1);
	}
	else {
		if (fp != NULL) {
			// ������ ���������
			fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

			// ������ ��� ����������
			fprintf(fp, "VARIABLES = x, y, z, Vx, Vy, Vz, Mag, Pressure , Normal , PAm, Zero\n");

			// ������ ���������� � �����
#if doubleintprecision == 1
			//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
			//if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
			fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
			//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
			//if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
			fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif
			

			integer i = 0; // �������� 
			//integer j = 0; // ����� for

			// ������ x
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ y
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ z
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");


			// ������ �������������� Vx ��������
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[VELOCITY_X_COMPONENT][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ ������������ Vy ��������
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[VELOCITY_Y_COMPONENT][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			fprintf(fp, "\n");

			// ������ Vz ��������
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[VELOCITY_Z_COMPONENT][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			fprintf(fp, "\n");


			// ������ ������ ��������
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", sqrt(potent[VELOCITY_X_COMPONENT][i] * potent[VELOCITY_X_COMPONENT][i] + potent[VELOCITY_Y_COMPONENT][i] * potent[VELOCITY_Y_COMPONENT][i] + potent[VELOCITY_Z_COMPONENT][i] * potent[VELOCITY_Z_COMPONENT][i]));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ ��������
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[PRESS][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ Normal
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", rhie_chow[0][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ PAm
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", potent[PAM][i]); // rhie_chow[1][i]
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ Rhi-Chow
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", rhie_chow[2][i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// ������ ���������� � ���������� �����
			for (i = 0; i < ncell; i++) {
#if doubleintprecision == 1
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			    }
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {

					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
				}
#else
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
				}
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {

					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
				}
#endif
				
			}

			
			fclose(fp); // �������� �����
			
			
			printf("file is successfully written...OK.\n");
			// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
			   //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL
		}
		else {
			printf("exporttecplotxy360_3D problem close file .\n");
			system("pause");
			exit(1);
		}
	}
} // exporttecplotxy360_3D

// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
void exporttecplotxy360T_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal* potent)
{
	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("ALICEFLOW0_03.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_03.PLT", "wb");
#endif
	
	// �������� ����� ��� ������.
	if ((err) != 0) {
		printf("Create File Error\n");
	}
	else {
		if (fp!=NULL) {

		// ������ ���������
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

		// ������ ��� ����������
		fprintf(fp, "VARIABLES = x, y, z, Temp\n");

#if doubleintprecision == 1
		// ������ ���������� � �����
		fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
		// ������ ���������� � �����
		fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif

		

		integer i = 0; // �������� 
		//integer j = 0; // ����� for

		// ������ x
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ y
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ z
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z));
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");


		// ������ �����������
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", potent[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ ���������� � ���������� �����
		for (i = 0; i < ncell; i++) {
#if doubleintprecision == 1
			//fprintf(fp, "%lld %lld %lld %lld ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
			//fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			fprintf(fp, "%d %d %d %d %d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i], nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#else
			//fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
			//fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			fprintf(fp, "%d %d %d %d %d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i], nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#endif

				}

		fclose(fp); // �������� �����
		printf("file is successfully written...OK.\n");
		// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
		   //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);
	}
	}
} // exporttecplotxy360T_3D

// ���������� ������ ����� ��������� ���������� 19N ����������� ������.
// ����� ����� ����������� ���������� ������������ ��������� �����.

// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 1 � 3.
void exporttecplotxy360T_3D_part1and3(TEMPER &t, integer maxelm, integer maxbound, bool bextendedprint, integer ncell,
									  int** nvtx, int** nvtxcell, TOCHKA* pa,
									  BOUND* border_neighbor, integer ivarexport, int** ptr_out)
{

	if (bvery_big_memory) {

		t.database.maxelm = maxelm;
		t.database.ncell = ncell;

		// extended print integer �� �������������.

		// ���� ������ ��� ���������� ����� �� � ���� ����������.
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

		// ���������� ����� ������� ������� � ������� ��������.
		t.database.x = (float*)malloc(maxelm*sizeof(float));
		t.database.y = (float*)malloc(maxelm*sizeof(float));
		t.database.z = (float*)malloc(maxelm*sizeof(float));
		t.database.ptr = new int*[2];
		for (integer j = 0; j < 2; j++) {
			t.database.ptr[j] = new int[maxelm];
		}
		for (integer j = 0; j < 2; j++) {
			for (integer i = 0; i < maxelm; i++) {
				t.database.ptr[j][i] = ptr_out[j][i];
			}
		}
		t.database.nvtxcell = new int*[8];
		for (integer i = 0; i < 8; i++) {
			t.database.nvtxcell[i] = new int[ncell];
		}

		// ������ x
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
		// ������ y
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

		// ������ z
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

		// ������ ���������� � ���������� �����
		for (integer i = 0; i < ncell; i++) {
			for (integer j = 0; j < 8; j++) {
				t.database.nvtxcell[j][i] = nvtxcell[j][i];
			}
		}

	}
	else {

		t.database.maxelm = 0;
		t.database.ncell = 0;
		//������ � ���� ����� ���������.
		t.database.x = nullptr;
		t.database.y = nullptr;
		t.database.z = nullptr;
		t.database.nvtxcell = nullptr;
		t.database.ptr = nullptr;


		// ����������� ������
		// ��� ����������� ������ �� �������� ����� � ��������� ����.
		// bextendedprint=true; ����������� ������. 

		// ivarexport == 1 ���������� ������ ���� ����������,
		// ivarexport == 2 ���������� ������ �������������,
		// ivarexport == 3 ���������� � ���� ���������� � �������������.

		FILE *fp=NULL;
		
#ifdef MINGW_COMPILLER
		int err = 0;
		fp = fopen64("ALICEFLOW0_06_temp_part1.txt", "w");
		if (fp == NULL) err = 1;
#else
		errno_t err = 0;
		err = fopen_s(&fp, "ALICEFLOW0_06_temp_part1.txt", "w");
#endif
		
		// �������� ����� ��� ������:
		// ���� ������� �� ��� ������: 
		// 1 � 3 ����� ������������ �����
		// ������ ����� � ������������ ������� ������������
		// ����� �������. ����� ���������� ������ ����� ������� � �����
		// ���������� ������ ������������ ����������� ������.
		// �������� ������ 19N.

		// ������ ������ 1 � 3
		if ((err) != 0) {
			printf("Create File temp part1 Error\n");
			//getchar();
			system("pause");

		}
		else {
			// ������ ���������
			/*fprintf(fp, "TITLE = \"ALICEFLOW0_06\"\n");

			// ������ ��� ����������
			switch (ivarexport) {
			case 1: fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n"); break; // ���������� ������ ���� ����������
			case 2: fprintf(fp, "VARIABLES = x, y, z, Speed, Pressure, Vx, Vy, Vz, Rho, Mu, Mut\n"); break; // ���������� ������ �������������
			case 3: fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, Speed, Pressure, Vx, Vy, Vz, Rho, Mu, Mut\n"); break; // ���������� � ����������� � �������������
			default: printf("Error export tecplot. Nonselected exporting variables...\n"); getchar();
			}

			#if doubleintprecision == 1
			    // ������ ���������� � �����
			    fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
			#else
			    // ������ ���������� � �����
			    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
			#endif

			
			*/

			integer i = 0; // �������� 
			//integer j=0; // ����� for

			// ������ x
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (border_neighbor[i].Norm) {// ��������� ���������� ������� � �������
					case E_SIDE: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI] - 1].x);
						break;
					case W_SIDE: fprintf(fp, "%+.16f ", pa[nvtx[1][border_neighbor[i].iI] - 1].x);
						break;
					case N_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					case S_SIDE:fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					case T_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					case B_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// ������ y
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (border_neighbor[i].Norm) {// ��������� ���������� ������� � �������
					case E_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					case W_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					case N_SIDE: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI] - 1].y);
						break;
					case S_SIDE: fprintf(fp, "%+.16f ", pa[nvtx[2][border_neighbor[i].iI] - 1].y);
						break;
					case T_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					case B_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// ������ z
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (border_neighbor[i].Norm) {// ��������� ���������� ������� � �������
					case E_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case W_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case N_SIDE: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case S_SIDE:fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case T_SIDE: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI] - 1].z);
						break;
					case B_SIDE: fprintf(fp, "%+.16f ", pa[nvtx[4][border_neighbor[i].iI] - 1].z);
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			fclose(fp); // �������� �����
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
				// ������ ���������� � ���������� �����
				for (i = 0; i < ncell; i++) {
#if doubleintprecision == 1
					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#else
					fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
#endif
					
				}

				fclose(fp); // �������� �����
				printf("export tecplot temperature part3 is successfully written...OK.\n");
			}
		}
	}

} // exporttecplotxy360T_3D_part1and3

// 10 ������ 2016 . �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
// � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2.
void exporttecplotxy360T_3D_part2binary(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint)
{
    // ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage=false;

	FILE *fp=NULL;
    FILE *fp1=NULL; // ����� 1 ��� 3
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp=fopen64("ALICEFLOW0_07_temp.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");
#endif
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	
	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
	if ((err) != 0) {
		printf("Create File temp Error in exporttecplotxy360T_3D_part2binary in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {
        
        int c; // �������� ������
		integer ivarexport=1; // �� ��������� ������ ���� ����������:
		integer i=0; // ������� �����

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

        
			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count>0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (i=0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport=3; // ������� ��� ����������� ������ �������
				}
			}
            
			if (ivarexport==1) {
				// ������ ���������
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


		         // ������ ��� ����������
		        // fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
				 char variablestit[] = "VARIABLES = x, y, z, Temp, Lam";
				 //31
				 for (i = 0; i <= 29; i++) {
					 symbol = variablestit[i];
					 fwrite(&symbol, sizeof(char), 1, fp);
				 }
				  symbol = '\n';
				 fwrite(&symbol, sizeof(char), 1, fp);


				 // ������ ���������� � �����
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
					 // extended printeger �� �������������.

					 // ������ x
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
					 // ������ y
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
					 // ������ z
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
				// ������ ���������
		        fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// ������ ����� ������� ������� � ������������� � �������������:
				fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z\n");
					
#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#endif

			
				if (bvery_big_memory) {
					// ������ x
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
					// ������ y
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
					// ������ z
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
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// ������ ������ �����
        
        // ������ ���� ���������� ������������ ������.
		

		// ������ �����������
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport==3) {
			char symbol;
			doublereal fnumber;
			// Speed
            for (i=0;i<maxelm; i++) {
				if (t.ptr[1][i]>-1) {
					doublereal svx=f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					doublereal svy=f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					doublereal svz=f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
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
				    	doublereal svx=f[idfluid].potent[VELOCITY_X_COMPONENT][i+maxelm]*f[idfluid].potent[VELOCITY_X_COMPONENT][i+maxelm];
					    doublereal svy=f[idfluid].potent[VELOCITY_Y_COMPONENT][i+maxelm]*f[idfluid].potent[VELOCITY_Y_COMPONENT][i+maxelm];
					    doublereal svz=f[idfluid].potent[VELOCITY_Z_COMPONENT][i+maxelm]*f[idfluid].potent[VELOCITY_Z_COMPONENT][i+maxelm];
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
                 //  fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
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
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i+maxelm]); // VX
					   fnumber = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
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
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
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
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i+maxelm]); // VY
					   fnumber = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
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
                  // fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
				   fnumber = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
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
                      // fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i+maxelm]); // VZ
					   fnumber = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
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
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					//--->//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
					fnumber = f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]];
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
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
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
					fnumber = f[0].prop_b[MU_DYNAMIC_VISCOSITY][i];
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

			// Mut // ������������ ������������ ��������
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

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if ((f[0].iflowregime== VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[0].iflowregime== VISCOSITY_MODEL::SMAGORINSKY)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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

			

			// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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

			// ������������� ����.
			for (i=0; i<t.maxelm+t.maxbound; i++) {
				Tx[i]=0.0;
				Ty[i]=0.0;
				Tz[i]=0.0;
			}

			// ���������� ����������.
			for (i=0; i<t.maxelm; i++) {
				// ������ ���������� ����.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, false, 
					t.border_neighbor,  Tx, Ty, Tz, t.ilevel_alice);
			}

			for (i=0; i<t.maxelm; i++) {
				// ������ ��������� ����.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, true, 
					t.border_neighbor,Tx, Ty, Tz, t.ilevel_alice);
			}
			

			// ���������� � ����.

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

			// ������������ ����������� ������.
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
			// ������ ���������� � ���������� �����
			for (i = 0; i < t.database.ncell; i++) {
#if doubleintprecision == 1
				//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
				//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);

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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}
		*/

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}

// 30.09.2019 ������ ��������� ����� � ������������ �������� � ������ ������������� �������.
void save_velocity_for_init(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count) {


	
	FILE* fp_inicialization_data=NULL;
#ifdef MINGW_COMPILLER
	int err_inicialization_data = 0;
	fp_inicialization_data=fopen64("load.txt", "w");
	if (fp_inicialization_data==NULL) err_inicialization_data = 1;
#else
	errno_t err_inicialization_data = 0;
	err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "w");
#endif

	
	if (err_inicialization_data != 0) {
		// �������� �������� ��� ���� �����������.
		printf("Create File load.txt Error\n");
		//getchar();
		system("pause");
	}
	else {

		if (b_on_adaptive_local_refinement_mesh) {
			// ������������� ���� ��������� �� ���� �����.


			// �������� ������ ��� ����� FLUID ����.
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
			int ifr = 0;
			switch (f[0].iflowregime) {
			case VISCOSITY_MODEL::LAMINAR: ifr = 0;
				break;
			case VISCOSITY_MODEL::ZEROEQMOD: ifr = 1;
				break;
			case VISCOSITY_MODEL::SMAGORINSKY: ifr = 2;
				break;
			case VISCOSITY_MODEL::RNG_LES: ifr = 3;
				break;
			case VISCOSITY_MODEL::RANS_SPALART_ALLMARES:
				ifr = 4;
				break;
			case VISCOSITY_MODEL::RANS_MENTER_SST:
				ifr = 5;
				break;
			case VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST:
				ifr = 5;
				break;
			case VISCOSITY_MODEL::RANS_STANDART_K_EPS:
				ifr = 6;
				break;
			default:
				ifr = 0;
				break;
			}


#if doubleintprecision == 1
			fprintf_s(fp_inicialization_data, "%lld\n %lld\n %lld\n", f[0].maxnod, f[0].maxelm, (integer)(ifr));
#else
			fprintf_s(fp_inicialization_data, "%d\n %d\n %d\n", f[0].maxnod, f[0].maxelm, ifr);
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



			// ��������� ���������� �������� �� ������� ����� � ���� �����.
			// ������� � �������� ������������.
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

				doublereal Speed = sqrt((f[0].potent[VELOCITY_X_COMPONENT][i])*(f[0].potent[VELOCITY_X_COMPONENT][i])+(f[0].potent[VELOCITY_Y_COMPONENT][i])*(f[0].potent[VELOCITY_Y_COMPONENT][i])+(f[0].potent[VELOCITY_Z_COMPONENT][i])*(f[0].potent[VELOCITY_Z_COMPONENT][i]));
				if (Speed > SpeedMax) SpeedMax = Speed;
				if (Speed < SpeedMin) SpeedMin = Speed;

				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
				volume3D(i, f[0].nvtx, f[0].pa, dx, dy, dz);
				
				if (Ux != nullptr) {
					Ux[inode1] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode2] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode3] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode4] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode5] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode6] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode7] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
					Ux[inode8] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_X_COMPONENT][i];
				}

				if (Uy != nullptr) {
					Uy[inode1] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode2] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode3] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode4] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode5] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode6] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode7] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
					Uy[inode8] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Y_COMPONENT][i];
				}

				if (Uz != nullptr) {
					Uz[inode1] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode2] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode3] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode4] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode5] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode6] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode7] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
					Uz[inode8] += 0.125*dx*dy*dz*f[0].potent[VELOCITY_Z_COMPONENT][i];
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
				Ux[inode1] += (1.0/(0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode2] += (1.0/(0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode3] += (1.0/(0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode4] += (1.0 / (0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode5] += (1.0 / (0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode6] += (1.0 / (0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode7] += (1.0 / (0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];
				Ux[inode8] += (1.0/(0.5*dx))*f[0].potent[VELOCITY_X_COMPONENT][i];

				vesaX[inode1] += (1.0 / (0.5*dx));
				vesaX[inode2] += (1.0 / (0.5*dx));
				vesaX[inode3] += (1.0 / (0.5*dx));
				vesaX[inode4] += (1.0 / (0.5*dx));
				vesaX[inode5] += (1.0 / (0.5*dx));
				vesaX[inode6] += (1.0 / (0.5*dx));
				vesaX[inode7] += (1.0 / (0.5*dx));
				vesaX[inode8] += (1.0 / (0.5*dx));

				Uy[inode1] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode2] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode3] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode4] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode5] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode6] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode7] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];
				Uy[inode8] += (1.0 / (0.5*dy))*f[0].potent[VELOCITY_Y_COMPONENT][i];

				vesaY[inode1] += (1.0 / (0.5*dy));
				vesaY[inode2] += (1.0 / (0.5*dy));
				vesaY[inode3] += (1.0 / (0.5*dy));
				vesaY[inode4] += (1.0 / (0.5*dy));
				vesaY[inode5] += (1.0 / (0.5*dy));
				vesaY[inode6] += (1.0 / (0.5*dy));
				vesaY[inode7] += (1.0 / (0.5*dy));
				vesaY[inode8] += (1.0 / (0.5*dy));


				Uz[inode1] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode2] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode3] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode4] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode5] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode6] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode7] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];
				Uz[inode8] += (1.0 / (0.5*dz))*f[0].potent[VELOCITY_Z_COMPONENT][i];

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

			// ���������������� �� ������ �������� ���������� ��������.
			for (integer i = 0; i < f[0].maxnod; i++) {
				if (fabs(vol[i]) < 1.0e-36) {
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

			//���������� �� �������������� ����������� ���������� ���� ���������� ����� 
			// � ���������. 21.03.2019
			Speed_coef /= SpeedMax;// ����������� ������ 1. (�������� 1.6).
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

			// ������������ ������������ ��������.
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

			// ������������ ������������ ��������.
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

			// ������� �� ���� ����� ��������.

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
				fprintf(fp_inicialization_data, "%e ", f[0].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
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
				fprintf(fp_inicialization_data, "%e ", f[0].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
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
				fprintf(fp_inicialization_data, "%e ", f[0].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
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
			int ifr = 0;
			switch (f[0].iflowregime) {
			case VISCOSITY_MODEL::LAMINAR: ifr = 0;
				break;
			case VISCOSITY_MODEL::ZEROEQMOD: ifr = 1;
				break;
			case VISCOSITY_MODEL::SMAGORINSKY: ifr = 2;
				break;
			case VISCOSITY_MODEL::RNG_LES: ifr = 3;
				break;
			case VISCOSITY_MODEL::RANS_SPALART_ALLMARES:
				ifr = 4;
				break;
			case VISCOSITY_MODEL::RANS_MENTER_SST:
				ifr = 5;
				break;
			case VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST:
				ifr = 5;
				break;
			case VISCOSITY_MODEL::RANS_STANDART_K_EPS:
				ifr = 6;
				break;
			default:
				ifr = 0;
				break;
			}

#if doubleintprecision == 1
		fprintf_s(fp_inicialization_data, "%lld\n %lld\n %lld\n", maxelm, ncell, (integer)(ifr));
#else
		fprintf_s(fp_inicialization_data, "%d\n %d\n %d\n", maxelm, ncell, ifr);
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
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
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
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
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
					fprintf_s(fp_inicialization_data, "%e ", f[0].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
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

// 10 ������ 2016 . �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
// � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2.
void exporttecplotxy360T_3D_part2_apparat_hot( integer maxelm, integer ncell,
	FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate,
	bool bextendedprint, integer ikey, BLOCK* &b)
{
	const bool lite_export = true;
	// 16 ������ ����� ������� ��������� ������ �������,
	// ������ ���������� ����� ������.

	printf("ionly_solid_visible =%d\n", ionly_solid_visible);	

	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // ����� 1 ��� 3
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	doublereal* temp_shadow = nullptr;
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
		total_deformation_shadow = new doublereal*[4];
		for (integer j_1 = 0; j_1 < 4; j_1++) {
			total_deformation_shadow[j_1] = new doublereal[t.maxelm + t.maxbound];
			for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
				total_deformation_shadow[j_1][i_1] = t.total_deformation[j_1][i_1];
			}
		}
	}


	static int inum_exp = -1;
	inum_exp++;

	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	int err = 0;
	// ���� ��������� �� �� �������� ������ �� �����.
	if (inum_exp>=0) {
	switch (ikey) {
		case 0:  fp=fopen64("ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
		case 1:  fp=fopen64("ALICEFLOW0_27_temp_apparat_hot.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp=fopen64("ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
    }
	}
	/*
	if (inum_exp == 1) {
		switch (ikey) {
		case 0:  fp = fopen64("ALICEFLOW0_07_temp_apparat_hot1.PLT", "wb");  break;
		case 1:  fp = fopen64("ALICEFLOW0_27_temp_apparat_hot1.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_apparat_hot1.PLT", "wb");  break;
		}
	}
	if (inum_exp == 2) {
		switch (ikey) {
		case 0:  fp = fopen64("ALICEFLOW0_07_temp_apparat_hot2.PLT", "wb");  break;
		case 1:  fp = fopen64("ALICEFLOW0_27_temp_apparat_hot2.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_apparat_hot2.PLT", "wb");  break;
		}
	}
	if (inum_exp == 3) {
		switch (ikey) {
		case 0:  fp = fopen64("ALICEFLOW0_07_temp_apparat_hot3.PLT", "wb");  break;
		case 1:  fp = fopen64("ALICEFLOW0_27_temp_apparat_hot3.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_apparat_hot3.PLT", "wb");  break;
		}
    }
	if (inum_exp == 4) {
		switch (ikey) {
		case 0:  fp = fopen64("ALICEFLOW0_07_temp_apparat_hot4.PLT", "wb");  break;
		case 1:  fp = fopen64("ALICEFLOW0_27_temp_apparat_hot4.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_apparat_hot4.PLT", "wb");  break;
		}
	}
	if (inum_exp >= 5) {
		switch (ikey) {
		case 0:  fp = fopen64("ALICEFLOW0_07_temp_apparat_hot5.PLT", "wb");  break;
		case 1:  fp = fopen64("ALICEFLOW0_27_temp_apparat_hot5.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_apparat_hot5.PLT", "wb");  break;
		}
	}
	*/
	if (fp == NULL) err = 1;
#else
	errno_t err;
	// ���� ��������� �� �� ������ �� ����� ����� ���������.
	if (inum_exp >= 0) {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_apparat_hot.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot.PLT", "wb");  break;
		}
	}
	/*
	if (inum_exp == 1) {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot1.PLT", "wb");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_apparat_hot1.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot1.PLT", "wb");  break;
		}
	}
	if (inum_exp >= 2) {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot2.PLT", "wb");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_apparat_hot2.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_apparat_hot2.PLT", "wb");  break;
		}
	}
	*/
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in exporttecplotxy360T_3D_part2_apparat_hot in my_export_tecplot.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // �������� ������
		integer ivarexport = 1; // �� ��������� ������ ���� ����������:
		integer i = 0; // ������� �����

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


			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count>0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // ������� ��� ����������� ������ �������
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
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// ������ ��� ����������
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
				if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
					//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux\n");  // ���������� ������ ���� ����������
					fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������

				}
				else {

					fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // ���������� ������ ���� ����������
				}
#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif

				

				if (bvery_big_memory) {
					// extended printeger �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
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
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ������ �� ������
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
					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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
							
							integer inode2W= t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
							integer inode3W= t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
							integer inode6W= t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
							integer inode7W= t.neighbors_for_the_internal_node[W_SIDE][0][inode8];
							
							
							integer inode5B= t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
							integer inode6B= t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
							integer inode7B= t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
							integer inode8B= t.neighbors_for_the_internal_node[B_SIDE][0][inode4];
							

							
							integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2]; 
							integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1]; 
							integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6]; 
							integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5]; 
							
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


							// ������������ ������ ������� ����.
							// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
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
									if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[j_4][inode5] = t.total_deformation[j_4][inode1];
											total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode2];
											total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode3];
											total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode4];
										}
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
									if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
										// total_deformation
										for (integer j_4 = 0; j_4 < 4; j_4++) {
											total_deformation_shadow[j_4][inode2] = t.total_deformation[j_4][inode1];
											total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode4];
											total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode5];
											total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode8];
										}
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
									if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
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
						if (ncell_shadow == 0) {
							// � ��� ������ ���� ������� ����, ������� �� ���������� ������  ��������.
							ncell_shadow = ncell;
							ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;
						}
					}
				}
				// ������ ����� ������� ������� � ������������� � �������������:
				//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
				fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				

#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif

				
				

				if (bvery_big_memory) {
					if (lite_export) {
						// ������ x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					else {
						// ������ x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ z
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
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// ������ ������ �����

		// ������ ���� ���������� ������������ ������.


		// ������ �����������
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%+.3f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.3f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%+.16f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.16f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
					doublereal svy = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
					doublereal svz = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // ������������ ������������ ��������
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]/ f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);
					
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
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]/ f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if ((f[0].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[0].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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

		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {

			doublereal *Tx = nullptr;
			doublereal *Ty = nullptr;
			doublereal *Tz = nullptr;
			Tx = new doublereal[t.maxelm + t.maxbound];
			Ty = new doublereal[t.maxelm + t.maxbound];
			Tz = new doublereal[t.maxelm + t.maxbound];

			// ������������� ����.
			for (i = 0; i<t.maxelm + t.maxbound; i++) {
				Tx[i] = 0.0;
				Ty[i] = 0.0;
				Tz[i] = 0.0;
			}

			// ���������� ����������.
			for (i = 0; i<t.maxelm; i++) {
				// ������ ���������� ����.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, false,
					t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
			}

			for (i = 0; i<t.maxelm; i++) {
				// ������ ��������� ����.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, true,
					t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
			}


			// ���������� � ����.

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
			

			// ������������ ����������� ������.
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

			}

			if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
				for (integer j_6 = 0; j_6 < 4; j_6++) {

					// ������ ������ ����������
					if (lite_export) {
						for (i = 0; i < maxelm; i++) {


							if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
								if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
							if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
								fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
							}
							else {
								fprintf(fp, "%e ", t.total_deformation[j_6][i]);
							}
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						if (bextendedprint) {
							for (i = maxelm; i < maxelm + t.maxbound; i++) {
								if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
			}

		if (bvery_big_memory) {
			// ������ ���������� � ���������� �����
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
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						
					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

						integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
						integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
						integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
						integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


						integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
						integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
						integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
						integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



						integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
						integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
						integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
						integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];

						
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

						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{ 
							
							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
							// ������������ ������� ���� � ��������.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// ������������ ������� ���� � ��������.
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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
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
// ����������������� �� ���� �����.
// ������ ������� ���������� ��� ���������� �������� ��� ������� �� ������ �������� �������.
// �.�. � ������ �� ������ ���������� �� ����������� � �������� ������� �� ������� ������������ �������������� ������������ �������� ����������� ����������� �� ����� ������� ��������.
void tecplot360patcher_for_print_in_report() {
	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) {
		// ������ ������������� ������.
		fp = fopen64("ALICEFLOW0_08_temp.PLT", "r");
    }
	else {
		fp = fopen64("ALICEFLOW0_07_temp.PLT", "r");
	}
	if (fp==NULL) err = 1;
#else
	errno_t err = 0;
	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) {
		// ������ ������������� ������.
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
		
#ifdef MINGW_COMPILLER
		int err1 = 0;
		fp1 = fopen64("ALICEFlow0_07_Visualisation_Magement.PLT", "w");
		if (fp1==NULL) err1 = 1;
#else
		errno_t err1 = 0;
		err1 = fopen_s(&fp1, "ALICEFlow0_07_Visualisation_Magement.PLT", "w");
#endif
		
		if ((err1) != 0) {
			printf("Create File ALICEFlow0_07_Visualisation_Magement.PLT Error\n");
			//getchar();
			system("pause");

		}
		else {

			integer iN = -1, iE = -1, iVar = -1;
			char str[600] = "\0";
			const unsigned int VARIABLE_MAX_COUNT = 100;
			char** name = new char*[VARIABLE_MAX_COUNT];
			for (integer j32 = 0; j32 < VARIABLE_MAX_COUNT; j32++) {
				name[j32] = new char[600]; // ����� �������.
			}

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
				for (size_t i = 0; i < length_str; i++) if (str[i] == ',') {
					if (iVar < VARIABLE_MAX_COUNT) {
						size_t j32 = i - 1;
						while (str[j32] != ' ') j32--;
						j32++;
						size_t j33 = 0;
						while (str[j32] != ',') {
							name[iVar][j33] = str[j32];
							j33++;
							j32++;
						}
						name[iVar][j33] = '\0';
					}
					else {
						printf("ERROR!!! Number of VARIABLES > 100.\n");
						system("pause");
					}
					iVar++;
					//printf("i=%d iVar=%lld\n",i, iVar);
					//system("pause");
				}
				if (iVar < VARIABLE_MAX_COUNT) {
					size_t  j32 = length_str - 1;
					if (str[j32] == ' ') {
						// ����������� �������.
						while (str[j32] == ' ') j32--;
					}
					while (str[j32] != ' ') j32--;
					j32++;
					size_t j33 = 0;
					while ((j32 < length_str) && (str[j32] != ' ') && (str[j32] != '\n')) {
						name[iVar][j33] = str[j32];
						j33++;
						j32++;
					}
					name[iVar][j33] = '\0';
				}
				else {
					printf("ERROR!!! Number of VARIABLES > 100.\n");
					system("pause");
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
				for (size_t i = 0; i < length_str; i++) if (str[i] == '=') {

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
						for (size_t j = i + 1; j < length_str; j++) {
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
						for (size_t j = i + 1; j < length_str; j++) {
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

			// ���������� ���� �������.
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
				printf("%s compleate %lld is %lld\n", name[i], i + 1, iVar);
			}

			//printf("1 - X, 2 - Y, 3 - Z\n");
			//char ch2 = getchar();
			char ch2 = '2'; // �������������� ���� Y
			switch (pfpir.idir) {
			case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: ch2 = '1'; break;
			case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: ch2 = '2'; break;
			case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: ch2 = '3'; break;
			default:
				printf("Unknown directional in variable pfpir.idir=%d\n", pfpir.idir);
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

			if (1) {
#ifndef NO_OPENGL_GLFW

			pa_opengl = new TOCHKA[iN];
			pa_render = new TOCHKA[iN];
			n_render = iN;
			for (int i = 0; i < n_render; i++) {
				pa_opengl[i].x = func[0][i];
				pa_opengl[i].y = func[1][i];
				pa_opengl[i].z = func[2][i]; 

				pa_render[i].x = halfScreenWidth + scale_all * func[0][i];
				pa_render[i].y = halfScreenHeight + scale_all * func[1][i];
				pa_render[i].z = -abbys + scale_all * func[2][i];
			}


			doublereal dmin = 1.0e30;
			doublereal dmax = -1.0e30;

			//iCURENT_FUNC_openGL == 0;
			for (int i37 = 0; i37 < iN; i37++) {
				if (func[3][i37] > dmax) {
					dmax = func[3][i37];
				}
				if (func[3][i37] < dmin) {
					dmin = func[3][i37];
				}
			}

			minimum_val_for_render_pic = dmin;
			maximum_val_for_render_pic = dmax;

			int** nvtx_2 = new int* [8];
			for (int i_42 = 0; i_42 < 8; i_42++) {
				nvtx_2[i_42] = new int[iE];
				for (int i_43 = 0; i_43 < iE; i_43++) {
					nvtx_2[i_42][i_43] = (int)(nvtx_1[i_43][i_42]);
				}
			}

			iCURENT_FUNC_openGL = 0;
			potent_array_opengl = new doublereal * [iVar-3];
			for (int i_42 = 3; i_42 < iVar; i_42++) {
				potent_array_opengl[i_42 - 3] = new doublereal[iN];
				for (int i_43 = 0; i_43 < iN; i_43++) {
					potent_array_opengl[i_42 - 3][i_43] = func[i_42][i_43];
				}
			}
			iNUMBER_FUNC_openGL = iVar - 3;

			DrawZbufferColor( (int)(iE), (int)(iN), nvtx_2);
			delete[] pa_opengl;
			delete[] pa_render;
			n_render = -1;
			for (int i_42 = 0; i_42 < 8; i_42++) {
				delete[] nvtx_2[i_42];
			}
			delete[] nvtx_2;

			for (int i_42 = 3; i_42 < iVar; i_42++) {
				delete[] potent_array_opengl[i_42 - 3];
			}
			delete[] potent_array_opengl;
			potent_array_opengl = nullptr;

#endif

			}


			// ���������� ������� �� ������ ����������� ������.
			{
				// xyplot
				// ��� ������������ ����� xyplotT.txt.

				FILE *fpPLT = NULL;
				FILE* fpPLT1 = NULL;
				
#ifdef MINGW_COMPILLER
				int  err_PLT = 0;
				fpPLT = fopen64("xyplotT.txt", "w");
				if (fpPLT == NULL) err_PLT = 1;
#else
				errno_t err_PLT = 0;
				err_PLT = fopen_s(&fpPLT, "xyplotT.txt", "w");
#endif

#ifdef MINGW_COMPILLER
				int  err_PLT1 = 0;
				fpPLT1 = fopen64("xyplotT1.PLT", "w");
				if (fpPLT1 == NULL) err_PLT1 = 1;
#else
				errno_t err_PLT1 = 0;
				err_PLT1 = fopen_s(&fpPLT1, "xyplotT1.PLT", "w");
#endif


				if ((err_PLT) != 0) {
					printf("Create File xyplot Error\n");
					//getchar();
					system("pause");

				}
				else if (iE > 0) {
					if (fpPLT != NULL)
					{
						// �������� ! ��������� ������� ����� ����� ������� ����� ��������� ����� � ��������� ������� ������ ����� ����� ���������������.

						//doublereal x = 2.0245e-3, y = 0.2268e-3, z = 0.0125e-3; // ����� ����� ������� �������� �����
						doublereal x = 3e-3, y = -0.00005e-3, z = 0.57425e-3;
						x = Tochka_position_X0_for_XY_Plot;
						y = Tochka_position_Y0_for_XY_Plot;
						z = Tochka_position_Z0_for_XY_Plot;


						TOCHKA sp0, sp1; // ������� [sp0, sp1].
						doublereal xmin = 1.0e30, xmax = -1.0e30;
						doublereal ymin = 1.0e30, ymax = -1.0e30;
						doublereal zmin = 1.0e30, zmax = -1.0e30;
						for (integer iP = 0; iP < iE; iP++) {
							TOCHKA p;
							p.x = 0;
							p.y = 0;
							p.z = 0;
							integer ie1 = nvtx_1[iP][0] - 1;
							integer ie2 = nvtx_1[iP][1] - 1;
							integer ie3 = nvtx_1[iP][2] - 1;
							integer ie4 = nvtx_1[iP][3] - 1;
							integer ie5 = nvtx_1[iP][4] - 1;
							integer ie6 = nvtx_1[iP][5] - 1;
							integer ie7 = nvtx_1[iP][6] - 1;
							integer ie8 = nvtx_1[iP][7] - 1;
							p.x = 0.125*(func[0][ie1] + func[0][ie2] + func[0][ie3] + func[0][ie4] + func[0][ie5] + func[0][ie6] + func[0][ie7] + func[0][ie8]);
							p.y = 0.125*(func[1][ie1] + func[1][ie2] + func[1][ie3] + func[1][ie4] + func[1][ie5] + func[1][ie6] + func[1][ie7] + func[1][ie8]);
							p.z = 0.125*(func[2][ie1] + func[2][ie2] + func[2][ie3] + func[2][ie4] + func[2][ie5] + func[2][ie6] + func[2][ie7] + func[2][ie8]);

							if (p.x < xmin) xmin = p.x;
							if (p.x > xmax) xmax = p.x;
							if (p.y < ymin) ymin = p.y;
							if (p.y > ymax) ymax = p.y;
							if (p.z < zmin) zmin = p.z;
							if (p.z > zmax) zmax = p.z;
						}

						integer iplane = YZ_PLANE; // ��������� ���������������� �����.
						switch (idirectional_for_XY_Plot) {
						case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
							iplane = YZ_PLANE; // ��������� ���������������� �����.
							sp0.y = sp1.y = y;
							sp0.z = sp1.z = z;
							sp0.x = xmin;
							sp1.x = xmax;
							break;
						case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: //XZ
							iplane = XZ_PLANE; // ��������� ���������������� �����.
							sp0.x = sp1.x = x;
							sp0.z = sp1.z = z;
							sp0.y = ymin;
							sp1.y = ymax;
							break;
						case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
							iplane = XY_PLANE; // ��������� ���������������� �����.
							sp0.y = sp1.y = y;
							sp0.x = sp1.x = x;
							sp0.z = zmin;
							sp1.z = zmax;
							break;
						default:
							iplane = YZ_PLANE; // ��������� ���������������� �����.
							sp0.y = sp1.y = y;
							sp0.z = sp1.z = z;
							sp0.x = xmin;
							sp1.x = xmax;
							break;
						}

						doublereal dopusk = 1.0e30;
						for (integer iP = 0; iP < iE; iP++) {
							TOCHKA p;
							p.x = 0;
							p.y = 0;
							p.z = 0;
							// ���������� ��������� ������ ��.
							integer ie1 = nvtx_1[iP][0] - 1;
							integer ie2 = nvtx_1[iP][1] - 1;
							integer ie3 = nvtx_1[iP][2] - 1;
							integer ie4 = nvtx_1[iP][3] - 1;
							integer ie5 = nvtx_1[iP][4] - 1;
							integer ie6 = nvtx_1[iP][5] - 1;
							integer ie7 = nvtx_1[iP][6] - 1;
							integer ie8 = nvtx_1[iP][7] - 1;
							p.x = 0.125*(func[0][ie1] + func[0][ie2] + func[0][ie3] + func[0][ie4] + func[0][ie5] + func[0][ie6] + func[0][ie7] + func[0][ie8]);
							p.y = 0.125*(func[1][ie1] + func[1][ie2] + func[1][ie3] + func[1][ie4] + func[1][ie5] + func[1][ie6] + func[1][ie7] + func[1][ie8]);
							p.z = 0.125*(func[2][ie1] + func[2][ie2] + func[2][ie3] + func[2][ie4] + func[2][ie5] + func[2][ie6] + func[2][ie7] + func[2][ie8]);

							doublereal dx = fabs(func[0][ie2] - func[0][ie1]);
							doublereal dy = fabs(func[1][ie3] - func[1][ie1]);
							doublereal dz = fabs(func[2][ie5] - func[2][ie1]);
							//printf("dx=%e dy=%e dz=%e\n",dx,dy,dz); // �� ����� 20.05.2020.
							//getchar();

							if (distance_Point_to_Segment2(p, sp0, sp1) < kR_size_cell * (dx*dx+dy*dy+dz*dz)) {
								// ����� ����������� �������.
								if (kR_size_cell * (dx*dx + dy*dy + dz*dz) < dopusk) {
									dopusk = 0.49 * (dx*dx + dy*dy + dz*dz);
								}

							}
						}

						SORT_XY_PLOT_FLOW* xy_arr = new SORT_XY_PLOT_FLOW[iE];

						integer icounter = -1;

						for (integer iP = 0; iP < iE; iP++) {
							TOCHKA p;
							p.x = 0;
							p.y = 0;
							p.z = 0;
							// ���������� ��������� ������ ��.
							integer ie1 = nvtx_1[iP][0] - 1;
							integer ie2 = nvtx_1[iP][1] - 1;
							integer ie3 = nvtx_1[iP][2] - 1;
							integer ie4 = nvtx_1[iP][3] - 1;
							integer ie5 = nvtx_1[iP][4] - 1;
							integer ie6 = nvtx_1[iP][5] - 1;
							integer ie7 = nvtx_1[iP][6] - 1;
							integer ie8 = nvtx_1[iP][7] - 1;
							p.x = 0.125*(func[0][ie1] + func[0][ie2] + func[0][ie3] + func[0][ie4] + func[0][ie5] + func[0][ie6] + func[0][ie7] + func[0][ie8]);
							p.y = 0.125*(func[1][ie1] + func[1][ie2] + func[1][ie3] + func[1][ie4] + func[1][ie5] + func[1][ie6] + func[1][ie7] + func[1][ie8]);
							p.z = 0.125*(func[2][ie1] + func[2][ie2] + func[2][ie3] + func[2][ie4] + func[2][ie5] + func[2][ie6] + func[2][ie7] + func[2][ie8]);

							if (distance_Point_to_Segment2(p, sp0, sp1) < dopusk) {
								// ����� ����������� �������.

								icounter++;
								for (integer i32 = 3; i32 < iVar; i32++) {
									doublereal func_i32= 0.125*(func[i32][ie1] + func[i32][ie2] + func[i32][ie3] + func[i32][ie4] + func[i32][ie5] + func[i32][ie6] + func[i32][ie7] + func[i32][ie8]);
									switch (i32) {
									case 3: xy_arr[icounter].Vx = func_i32;  break;
									case 4: xy_arr[icounter].Vy = func_i32;  break;
									case 5: xy_arr[icounter].Vz = func_i32;  break;
									case 6: xy_arr[icounter].Speed = func_i32;  break;
									case 7: xy_arr[icounter].Pam = func_i32;  break;
									case 8: xy_arr[icounter].Press = func_i32;  break;
									case 9: xy_arr[icounter].Fbuf = func_i32;  break;
									case 10: xy_arr[icounter].GradPress = func_i32;  break;
									case 11: xy_arr[icounter].GradPam = func_i32;  break;
									case 12: xy_arr[icounter].massflux_in_gran = func_i32;  break;
									case 13: xy_arr[icounter].GradXVx = func_i32;  break;
									case 14: xy_arr[icounter].GradYVx = func_i32;  break;
									case 15: xy_arr[icounter].GradZVx = func_i32;  break;
									case 16: xy_arr[icounter].GradXVy = func_i32;  break;
									case 17: xy_arr[icounter].GradYVy = func_i32;  break;
									case 18: xy_arr[icounter].GradZVy = func_i32;  break;
									case 19: xy_arr[icounter].GradXVz = func_i32;  break;
									case 20: xy_arr[icounter].GradYVz = func_i32;  break;
									case 21: xy_arr[icounter].GradZVz = func_i32;  break;
									case 22: xy_arr[icounter].Curl = func_i32;  break;
									case 23: xy_arr[icounter].GradXPress = func_i32;  break;
									case 24: xy_arr[icounter].GradYPress = func_i32;  break;
									case 25: xy_arr[icounter].GradZPress = func_i32;  break;
									case 26: xy_arr[icounter].Nusha = func_i32;  break;
									case 27: xy_arr[icounter].GradXNusha = func_i32;  break;
									case 28: xy_arr[icounter].GradYNusha = func_i32;  break;
									case 29: xy_arr[icounter].GradZNusha = func_i32;  break;
									case 30: xy_arr[icounter].Ke = func_i32;  break;
									case 31: xy_arr[icounter].Omega = func_i32;  break;
									case 32: xy_arr[icounter].GradXKe = func_i32;  break;
									case 33: xy_arr[icounter].GradYKe = func_i32;  break;
									case 34: xy_arr[icounter].GradZKe = func_i32;  break;
									case 35: xy_arr[icounter].GradXOmega = func_i32;  break;
									case 36: xy_arr[icounter].GradYOmega = func_i32;  break;
									case 37: xy_arr[icounter].GradZOmega = func_i32;  break;
									case 38: xy_arr[icounter].F1 = func_i32; break;
									case 39: xy_arr[icounter].F2 = func_i32; break;
									case 40: xy_arr[icounter].F3 = func_i32; break;
									case 41: xy_arr[icounter].F4 = func_i32; break;
									case 42: xy_arr[icounter].F5 = func_i32; break;
									case 43: xy_arr[icounter].F6 = func_i32; break;
									case 44: xy_arr[icounter].F7 = func_i32; break;
									case 45: xy_arr[icounter].F8 = func_i32; break;
									case 46: xy_arr[icounter].F9 = func_i32; break;
									case 47: xy_arr[icounter].F10 = func_i32; break;
									case 48: xy_arr[icounter].F11 = func_i32; break;
									case 49: xy_arr[icounter].F12 = func_i32; break;
									case 50: xy_arr[icounter].F13 = func_i32; break;
									case 51: xy_arr[icounter].F14 = func_i32; break;
									case 52: xy_arr[icounter].F15 = func_i32; break;
									case 53: xy_arr[icounter].F16 = func_i32; break;
									case 54: xy_arr[icounter].F17 = func_i32; break;
									case 55: xy_arr[icounter].F18 = func_i32; break;
									case 56: xy_arr[icounter].F19 = func_i32; break;
									case 57: xy_arr[icounter].F20 = func_i32; break;
									case 58: xy_arr[icounter].F21 = func_i32; break;
									case 59: xy_arr[icounter].F22 = func_i32; break;
									case 60: xy_arr[icounter].F23 = func_i32; break;
									case 61: xy_arr[icounter].F24 = func_i32; break;
									case 62: xy_arr[icounter].F25 = func_i32; break;
									case 63: xy_arr[icounter].F26 = func_i32; break;
									case 64: xy_arr[icounter].F27 = func_i32; break;
									case 65: xy_arr[icounter].F28 = func_i32; break;
									case 66: xy_arr[icounter].F29 = func_i32; break;
									case 67: xy_arr[icounter].F30 = func_i32; break;
									case 68: xy_arr[icounter].F31 = func_i32; break;
									case 69: xy_arr[icounter].F32 = func_i32; break;
									default :
										printf("ERROR!!! number of VARIABLES > 69. XYPLOT diagnostic.\n");
										system("PAUSE");
										break;
									}
									 
								}
								

								switch (idirectional_for_XY_Plot) {
								case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
									xy_arr[icounter].x_argument = p.x;
									break;
								case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: // XZ
									xy_arr[icounter].x_argument = p.y;
									break;
								case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
									xy_arr[icounter].x_argument = p.z;
									break;
								default:
									// YZ
									xy_arr[icounter].x_argument = p.x;
									break;
								}

							}
						}

						// ���������� �� �����������.
						std::sort(xy_arr, xy_arr + icounter + 1, compare_XY_PLOT_FLOW);

						// ���� ��������� �������� � ����� � ��� �� 
						// ����������, �� �� ��������� ��������� �������,
						// �������� ������ � ����� ������ ��������.
						doublereal pos = xy_arr[0].x_argument;
						integer jposY = 0, jmaximumPos = jposY;
						bool bOkTemp = true;
						if ((iVar >= 3) && (name[3][0] == 'T') && (name[3][1] == 'e') && (name[3][2] == 'm') && (name[3][3] == 'p')) {
							// ��� ��� � ������������� ����� ��������  x, y, z �������� ����������� ������� ��������� � 
							// xy_arr[jposY].Vx;// �����������.
							bOkTemp = true;
						}
						else {
							printf("problem Temp not found in PLT file.");
							system("pause");
							bOkTemp = false; // ����������������� �� ������ ��������, ��� ��� ����������� ���.
						}
						doublereal maximumTemp = xy_arr[jposY].Vx;// �����������.
						if (0) {
							for (integer j = 1; j <= icounter; j++) {
								if (pos == xy_arr[j].x_argument) {
									xy_arr[j] = xy_arr[jposY];
								}
								else {
									jposY = j;
									pos = xy_arr[j].x_argument;
								}
							}
						}
						else {
							for (integer j = 1; j <= icounter; j++) {
								if (pos == xy_arr[j].x_argument) {
									if (xy_arr[j].Vx > maximumTemp) {
										maximumTemp = xy_arr[j].Vx;
										jmaximumPos = j;
									}
								}
								else {
									// ������ ����������� �������� �����������
									// � �����������, ��� ���� �������������.
									for (integer j31 = jposY; j31 < j; j31++) {
										xy_arr[j31] = xy_arr[jmaximumPos];
									}
									jposY = j;
									maximumTemp = xy_arr[jposY].Vx;// �����������.
									jmaximumPos = jposY;									
									pos = xy_arr[j].x_argument;
								}
							}
						}

						fprintf(fpPLT, "Xo=%e, Yo=%e, Zo=%e, ", x, y, z);
						switch (idirectional_for_XY_Plot) {
						case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
							fprintf(fpPLT, "directional=X\n");
							break;
						case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: //XZ
							fprintf(fpPLT, "directional=Y\n");
							break;
						case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
							fprintf(fpPLT, "directional=Z\n");
							break;
						default: // YZ 
							fprintf(fpPLT, "directional=X\n");
							break;
						}

						fprintf(fpPLT, "position,");
						for (integer i32 = 3; i32 < iVar; i32++) {
							fprintf(fpPLT, "\t%s,",name[i32]);
						}
						fprintf(fpPLT, "\n");

						fprintf(fpPLT1, "VARIABLES= position %s\n", name[3]);

						for (integer j = 0; j <= icounter; j++) {

							fprintf(fpPLT1, "%+.16f %+.16f\n", xy_arr[j].x_argument, xy_arr[j].Vx);

							fprintf(fpPLT, "%+.16f ", xy_arr[j].x_argument);
							for (integer i32 = 3; i32 < iVar; i32++) {
								switch (i32) {
								case 3: fprintf(fpPLT, "%+.16f ", xy_arr[j].Vx); break;
								case 4: fprintf(fpPLT, "%+.16f ", xy_arr[j].Vy);  break;
								case 5: fprintf(fpPLT, "%+.16f ", xy_arr[j].Vz);  break;
								case 6: fprintf(fpPLT, "%+.16f ", xy_arr[j].Speed);   break;
								case 7: fprintf(fpPLT, "%+.16f ", xy_arr[j].Pam);   break;
								case 8: fprintf(fpPLT, "%+.16f ", xy_arr[j].Press);   break;
								case 9: fprintf(fpPLT, "%+.16f ", xy_arr[j].Fbuf);   break;
								case 10: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradPress);   break;
								case 11: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradPam);   break;
								case 12: fprintf(fpPLT, "%+.16f ", xy_arr[j].massflux_in_gran);   break;
								case 13: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradXVx);   break;
								case 14: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradYVx);   break;
								case 15: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradZVx);   break;
								case 16: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradXVy);   break;
								case 17: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradYVy);   break;
								case 18: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradZVy);   break;
								case 19: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradXVz);   break;
								case 20: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradYVz);   break;
								case 21: fprintf(fpPLT, "%+.16f ", xy_arr[j].GradZVz);   break;
								case 22: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].Curl);   break;
								case 23: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradXPress);   break;
								case 24: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradYPress);   break;
								case 25: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradZPress);   break;
								case 26: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].Nusha);  break;
								case 27: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradXNusha);   break;
								case 28: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradYNusha);   break;
								case 29: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradZNusha);   break;
								case 30: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].Ke);   break;
								case 31: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].Omega);   break;
								case 32: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradXKe);   break;
								case 33: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradYKe);  break;
								case 34: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradZKe);  break;
								case 35: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradXOmega);   break;
								case 36: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradYOmega);   break;
								case 37: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].GradZOmega);   break;
								case 38: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F1);   break;
								case 39: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F2);   break;
								case 40: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F3);   break;
								case 41: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F4);   break;
								case 42: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F5);   break;
								case 43: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F6);   break;
								case 44: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F7);   break;
								case 45: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F8);   break;
								case 46: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F9);   break;
								case 47: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F10);   break;
								case 48: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F11);   break;
								case 49: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F12);   break;
								case 50: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F13);   break;
								case 51: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F14);   break;
								case 52: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F15);   break;
								case 53: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F16);   break;
								case 54: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F17);   break;
								case 55: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F18);   break;
								case 56: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F19);   break;
								case 57: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F20);   break;
								case 58: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F21);   break;
								case 59: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F22);   break;
								case 60: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F23);   break;
								case 61: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F24);   break;
								case 62: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F25);   break;
								case 63: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F26);   break;
								case 64: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F27);   break;
								case 65: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F28);   break;
								case 66: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F29);   break;
								case 67: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F30);   break;
								case 68: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F31);   break;
								case 69: fprintf(fpPLT, "%+.16f ", xy_arr[icounter].F32);   break;
								default:
									printf("ERROR!!! number of VARIABLES > 69. XYPLOT diagnostic.\n");
									system("PAUSE");
									break;
								}
							}
							fprintf(fpPLT, "\n");
						}
						
						delete[] xy_arr;

						printf("XYPLOT graphics is construct successful...\n");

						fclose(fpPLT);
						fclose(fpPLT1);

					}
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
				printf("%s compleate %lld is %lld\n", name[i], i + 1, iVar);
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

			// ������������ ����������� ������.
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

			if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) {
				printf("file ALICEFLOW0_08_temp.PLT preobrazovan. Ok.\n");
			}
			else {
				printf("file ALICEFLOW0_07_temp.PLT preobrazovan. Ok.\n");
			}

			for (integer j32 = 0; j32 < VARIABLE_MAX_COUNT; j32++) {
				delete[] name[j32];
			}
			delete[] name;
			system("pause");

		}
	}

}// tecplot360patcher_for_print_in_report

// ��� ������� ����������� �� ����������.
void green_gauss_Stress(integer iP,
	doublereal**& potent, int**& nvtx, TOCHKA*& pa,
	int***& neighbors_for_the_internal_node, integer maxelm, bool bbond,
	BOUND*& border_neighbor, integer* ilevel_alice, integer iDATA, integer iTARGET, LINE_DIRECTIONAL iDIRECTIONAL);


void print_tetrahedron(FILE* &fp, bool b1bvisible, bool b2bvisible, bool b3bvisible, bool b4bvisible, 
	bool b5bvisible, bool b6bvisible, bool b7bvisible, bool b8bvisible,
	integer nvtxcell0, integer nvtxcell1, integer nvtxcell2, integer nvtxcell3, 
	integer nvtxcell4, integer nvtxcell5, integer nvtxcell6, integer nvtxcell7) {


	// ��������� ������ XOY plane

	integer buf = nvtxcell2;
	nvtxcell2 = nvtxcell3;
	nvtxcell3 = buf;

	buf = nvtxcell6;
	nvtxcell6 = nvtxcell7;
	nvtxcell7 = buf;

	bool bbuf = b3bvisible;
	b3bvisible = b4bvisible;
	b4bvisible = bbuf;

    bbuf = b7bvisible;
	b7bvisible = b8bvisible;
	b8bvisible = bbuf;


	if (((!b1bvisible)&&(!b5bvisible) && (b2bvisible) && (b3bvisible) && (b4bvisible) && (b6bvisible) && (b7bvisible) && (b8bvisible))||
		((!b3bvisible)&&(!b7bvisible) && (b1bvisible) && (b2bvisible) && (b4bvisible) && (b5bvisible) && (b6bvisible) && (b8bvisible))||
		((b1bvisible)&&(b2bvisible)&&(b3bvisible)&&(b4bvisible)&&(b5bvisible) && (b6bvisible) && (b7bvisible) && (b8bvisible)))
	{

		// ������ ������������ ����� ����������.
		if ((b1bvisible) && (b2bvisible) && (b3bvisible) && (b5bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell0, nvtxcell1, nvtxcell2, nvtxcell4);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell0, nvtxcell1, nvtxcell2, nvtxcell4);
#endif
		}

		if ((b2bvisible) && (b5bvisible) && (b6bvisible) && (b7bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell1, nvtxcell4, nvtxcell5, nvtxcell6);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell1, nvtxcell4, nvtxcell5, nvtxcell6);
#endif
		}

		if ((b2bvisible) && (b3bvisible) && (b7bvisible) && (b5bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell1, nvtxcell2, nvtxcell6, nvtxcell4);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell1, nvtxcell2, nvtxcell6, nvtxcell4);
#endif
		}

		if ((b3bvisible) && (b2bvisible) && (b4bvisible) && (b7bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell2, nvtxcell1, nvtxcell3, nvtxcell6);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell2, nvtxcell1, nvtxcell3, nvtxcell6);
#endif
		}

		if ((b7bvisible) && (b8bvisible) && (b6bvisible) && (b4bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell6, nvtxcell7, nvtxcell5, nvtxcell3);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell6, nvtxcell7, nvtxcell5, nvtxcell3);
#endif
		}

		if ((b7bvisible) && (b6bvisible) && (b2bvisible) && (b4bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell6, nvtxcell5, nvtxcell1, nvtxcell3);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell6, nvtxcell5, nvtxcell1, nvtxcell3);
#endif
		}
	}
	else if (((!b2bvisible) && (!b6bvisible) && (b1bvisible) && (b3bvisible) && (b4bvisible) && (b5bvisible) && (b7bvisible) && (b8bvisible)) || 
		((!b4bvisible) && (!b8bvisible) && (b1bvisible) && (b2bvisible) && (b3bvisible) && (b5bvisible) && (b6bvisible) && (b7bvisible)) ||
		((b1bvisible) && (b2bvisible) && (b3bvisible) && (b4bvisible) && (b5bvisible) && (b6bvisible) && (b7bvisible) && (b8bvisible)))
	{

		// ������ ������������ ����� ����������.
		if ((b2bvisible) && (b4bvisible) && (b1bvisible) && (b6bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell1, nvtxcell3, nvtxcell0, nvtxcell5);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell1, nvtxcell3, nvtxcell0, nvtxcell5);
#endif
		}

		if ((b4bvisible) && (b6bvisible) && (b8bvisible) && (b5bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell3, nvtxcell5, nvtxcell7, nvtxcell4);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell3, nvtxcell5, nvtxcell7, nvtxcell4);
#endif
		}

		if ((b4bvisible) && (b1bvisible) && (b5bvisible) && (b6bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell3, nvtxcell0, nvtxcell4, nvtxcell5);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell3, nvtxcell0, nvtxcell4, nvtxcell5);
#endif
		}

		if ((b1bvisible) && (b4bvisible) && (b3bvisible) && (b5bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell0, nvtxcell3, nvtxcell2, nvtxcell4);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell0, nvtxcell3, nvtxcell2, nvtxcell4);
#endif
		}

		if ((b5bvisible) && (b7bvisible) && (b8bvisible) && (b3bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell4, nvtxcell6, nvtxcell7, nvtxcell2);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell4, nvtxcell6, nvtxcell7, nvtxcell2);
#endif
		}

		if ((b5bvisible) && (b8bvisible) && (b4bvisible) && (b3bvisible)) {
#if doubleintprecision == 1
			fprintf(fp, "%lld %lld %lld %lld\n", nvtxcell4, nvtxcell7, nvtxcell3, nvtxcell2);
#else
			fprintf(fp, "%d %d %d %d\n", nvtxcell4, nvtxcell7, nvtxcell3, nvtxcell2);
#endif
		}

	}

}



// 10 ������ 2016 . �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
// � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2.
void exporttecplotxy360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t,
	integer flow_interior_count, integer ianimate, bool bextendedprint,
	integer ikey, BLOCK* &b, integer &lb)
{

	bool bTETRAHEDRON = false; // 10.09.2020

	

	bool bFLAGitoa = false;

	const bool lite_export = true;
	// 16 ������ ����� ������� ��������� ������ �������,
	// ������ ���������� ����� ������.

	printf("ionly_solid_visible =%d\n", ionly_solid_visible);

	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage = false;

	FILE *fp=NULL; // �������� ����.
	FILE *fp1=NULL; // ����� 1 ��� 3
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	doublereal* temp_shadow = nullptr;
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}


	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
	(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL)||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
	{

		// gamma_xy
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 2, 4, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 2, 4, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 1, 5, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 1, 5, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			
			t.total_deformation[STRAIN_XY][i_1] = t.total_deformation[4][i_1] + t.total_deformation[5][i_1];
			
		}

		// gamma_yz

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 3, 4, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 3, 4, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 2, 5, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 2, 5, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			t.total_deformation[STRAIN_YZ][i_1] = t.total_deformation[4][i_1] + t.total_deformation[5][i_1];

		}


		// gamma_zx
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 1, 4, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 1, 4, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 3, 5, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 3, 5, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			t.total_deformation[STRAIN_ZX][i_1] = t.total_deformation[4][i_1] + t.total_deformation[5][i_1];

		}

		// epsilon_x
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, XDEFORMATION, STRAIN_X, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, XDEFORMATION, STRAIN_X, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

		

		// epsilon_y
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, YDEFORMATION, STRAIN_Y, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, YDEFORMATION, STRAIN_Y, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

		// epsilon_z
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ���������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, ZDEFORMATION, STRAIN_Z, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// ������ ��������� ����.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, ZDEFORMATION, STRAIN_Z, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

		for (integer i = 0; i < t.maxelm+t.maxbound; i++) {
			doublereal beta_t_solid_x; 
			doublereal beta_t_solid_y; 
			doublereal beta_t_solid_z;

			if (i < t.maxelm) {
				beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][i] * t.prop[BETA_T_MECHANICAL][i];// ����������� ��������� ��������� ���������� 1/K.
				beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][i] * t.prop[BETA_T_MECHANICAL][i];
				beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][i] * t.prop[BETA_T_MECHANICAL][i];
			}
			else {
				// ��������� ��.
				integer i_b = i - t.maxelm;
				integer ii = t.border_neighbor[i_b].iI;
				beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][ii] * t.prop[BETA_T_MECHANICAL][ii];// ����������� ��������� ��������� ���������� 1/K.
				beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][ii] * t.prop[BETA_T_MECHANICAL][ii];
				beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][ii] * t.prop[BETA_T_MECHANICAL][ii];
			}


			//printf("%e %e %e %e %e\n", beta_t_solid_x, beta_t_solid_y, beta_t_solid_z, (t.potent[i] - t.operatingtemperature_copy), t.operatingtemperature_copy);
			//getchar();

			if (eqin.itemper > 0) {

				t.total_deformation[STRAIN_X][i] -= beta_t_solid_x * (t.potent[i] - t.operatingtemperature_copy);
		        t.total_deformation[STRAIN_Y][i] -= beta_t_solid_y * (t.potent[i] - t.operatingtemperature_copy);
				t.total_deformation[STRAIN_Z][i] -= beta_t_solid_z * (t.potent[i] - t.operatingtemperature_copy);
			}
		}

		//double **Dirichlet = new doublereal*[6];
		//for (integer i_11 = 0; i_11 < 6; i_11++) {
			//Dirichlet[i_11] = new doublereal[6];
		//}
		for (integer i_1 = 0; i_1 < t.maxelm+t.maxbound; i_1++) {


			doublereal E;
			doublereal nu;

			if (i_1 < t.maxelm) {
				E = t.prop[YOUNG_MODULE][i_1];
				nu = t.prop[POISSON_RATIO][i_1];
			}
			else {
				// ��������� ��.
				integer i_b = i_1 - t.maxelm;
				integer ii = t.border_neighbor[i_b].iI;

				E = t.prop[YOUNG_MODULE][ii];
				nu = t.prop[POISSON_RATIO][ii];
			}

			//doublereal beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1]; // ������������ ����, ����������� ��������� ��������� ����������.
			//doublereal beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][i_1] * t.prop[BETA_T_MECHANICAL][i_1];// ����������� ��������� ��������� ���������� 1/K.
			//doublereal beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
			//doublereal beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
			//doublereal Ex = t.prop[MULT_YOUNG_MODULE_X][i_1] * t.prop[YOUNG_MODULE][i_1]; // ������ ���� ��.
			//doublereal Ey = t.prop[MULT_YOUNG_MODULE_Y][i_1] * t.prop[YOUNG_MODULE][i_1]; // ������ ���� ��.
			//doublereal Ez = t.prop[MULT_YOUNG_MODULE_Z][i_1] * t.prop[YOUNG_MODULE][i_1]; // ������ ���� ��.
			//doublereal E = t.prop[YOUNG_MODULE][i_1];
			//doublereal nuyz = t.prop[MULT_POISSON_RATIO_YZ][i_1] * t.prop[POISSON_RATIO][i_1];
			//doublereal nuxz = t.prop[MULT_POISSON_RATIO_XZ][i_1] * t.prop[POISSON_RATIO][i_1];
			//doublereal nuxy = t.prop[MULT_POISSON_RATIO_XY][i_1] * t.prop[POISSON_RATIO][i_1];
			//doublereal nuzy = t.prop[MULT_POISSON_RATIO_ZY][i_1] * t.prop[POISSON_RATIO][i_1];
			//doublereal nuzx = t.prop[MULT_POISSON_RATIO_ZX][i_1] * t.prop[POISSON_RATIO][i_1];
			//doublereal nuyx = t.prop[MULT_POISSON_RATIO_YX][i_1] * t.prop[POISSON_RATIO][i_1];

			//doublereal nu = t.prop[POISSON_RATIO][i_1];
			
			/*
			doublereal Gxy, Gyz, Gxz;
			if (!t.bActiveShearModule[i_1]) {
				Gxy = Gyz = Gxz = Ex / (2.0 * (1.0 + nuxy));
			}
			else {
				Gyz = t.prop[SHEAR_MODULE_YZ][i_1];
				Gxz = t.prop[SHEAR_MODULE_XZ][i_1];
				Gxy = t.prop[SHEAR_MODULE_XY][i_1];
			}
			*/
			/*
			Dirichlet[0][0] = (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz
				* nuzx + nuyz * nuzy - 1.0)*Ex;
			Dirichlet[0][1] = -(nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ey;
			Dirichlet[0][2] = -(nuxy*nuyz + nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ez;
			Dirichlet[0][3] = 0.0;
			Dirichlet[0][4] = 0.0;
			Dirichlet[0][5] = 0.0;
			Dirichlet[1][0] = -(nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ex;
			Dirichlet[1][1] = (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz
				* nuzx + nuyz * nuzy - 1.0)*Ey;
			Dirichlet[1][2] = -(nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ez;
			Dirichlet[1][3] = 0.0;
			Dirichlet[1][4] = 0.0;
			Dirichlet[1][5] = 0.0;
			Dirichlet[2][0] = -(nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ex;
			Dirichlet[2][1] = -(nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ey;
			Dirichlet[2][2] = (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz
				* nuzx + nuyz * nuzy - 1.0)*Ez;
			Dirichlet[2][3] = 0.0;
			Dirichlet[2][4] = 0.0;
			Dirichlet[2][5] = 0.0;
			Dirichlet[3][0] = 0.0;
			Dirichlet[3][1] = 0.0;
			Dirichlet[3][2] = 0.0;
			Dirichlet[3][3] = Gxy;
			Dirichlet[3][4] = 0.0;
			Dirichlet[3][5] = 0.0;
			Dirichlet[4][0] = 0.0;
			Dirichlet[4][1] = 0.0;
			Dirichlet[4][2] = 0.0;
			Dirichlet[4][3] = 0.0;
			Dirichlet[4][4] = Gyz;
			Dirichlet[4][5] = 0.0;
			Dirichlet[5][0] = 0.0;
			Dirichlet[5][1] = 0.0;
			Dirichlet[5][2] = 0.0;
			Dirichlet[5][3] = 0.0;
			Dirichlet[5][4] = 0.0;
			Dirichlet[5][5] = Gxz;
			*/

			// Compute 3D constitutive matrix (linear continuum mechanics)
			/*
			doublereal Dirichlet[6][6] =
			{
				{1.0 - nu, nu, nu, 0.0, 0.0, 0.0},
				{nu, 1.0 - nu, nu, 0.0, 0.0, 0.0},
				{nu, nu, 1.0 - nu, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, (1.0 - 2.0 * nu) / 2.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 0.0, (1.0 - 2.0 * nu) / 2.0, 0.0},
				{0.0, 0.0, 0.0, 0.0, 0.0, (1.0 - 2.0 * nu) / 2.0}
			};

			for (int i_r = 0; i_r < 6; i_r++) {
				for (int i_2 = 0; i_2 < 6; i_2++) {
					Dirichlet[i_r][i_2] *= E / ((1.0 + nu) * (1.0 - 2.0 * nu));
				}
			}
			*/


			doublereal Gxy, Gyz, Gxz;
			doublereal Ex, Ey, Ez;
			doublereal nuyz, nuxz, nuxy, nuzy, nuzx, nuyx;

			if (i_1 < t.maxelm) {
				if (!t.bActiveShearModule[i_1]) {
					Gxy = Gyz = Gxz = E / (2.0 * (1.0 + nu));
				}
				else {
					Gyz = t.prop[SHEAR_MODULE_YZ][i_1];
					Gxz = t.prop[SHEAR_MODULE_XZ][i_1];
					Gxy = t.prop[SHEAR_MODULE_XY][i_1];
				}

				 Ex = t.prop[MULT_YOUNG_MODULE_X][i_1] * t.prop[YOUNG_MODULE][i_1]; // ������ ���� ��.
				 Ey = t.prop[MULT_YOUNG_MODULE_Y][i_1] * t.prop[YOUNG_MODULE][i_1]; // ������ ���� ��.
				 Ez = t.prop[MULT_YOUNG_MODULE_Z][i_1] * t.prop[YOUNG_MODULE][i_1]; // ������ ���� ��.

				 nuyz = t.prop[MULT_POISSON_RATIO_YZ][i_1] * t.prop[POISSON_RATIO][i_1];
				 nuxz = t.prop[MULT_POISSON_RATIO_XZ][i_1] * t.prop[POISSON_RATIO][i_1];
				 nuxy = t.prop[MULT_POISSON_RATIO_XY][i_1] * t.prop[POISSON_RATIO][i_1];
				 nuzy = t.prop[MULT_POISSON_RATIO_ZY][i_1] * t.prop[POISSON_RATIO][i_1];
				 nuzx = t.prop[MULT_POISSON_RATIO_ZX][i_1] * t.prop[POISSON_RATIO][i_1];
				 nuyx = t.prop[MULT_POISSON_RATIO_YX][i_1] * t.prop[POISSON_RATIO][i_1];
			}
			else {
				integer i_b = i_1 - t.maxelm;
				integer ii = t.border_neighbor[i_b].iI;

				if (!t.bActiveShearModule[ii]) {
					Gxy = Gyz = Gxz = E / (2.0 * (1.0 + nu));
				}
				else {
					Gyz = t.prop[SHEAR_MODULE_YZ][ii];
					Gxz = t.prop[SHEAR_MODULE_XZ][ii];
					Gxy = t.prop[SHEAR_MODULE_XY][ii];
				}

				Ex = t.prop[MULT_YOUNG_MODULE_X][ii] * t.prop[YOUNG_MODULE][ii]; // ������ ���� ��.
				Ey = t.prop[MULT_YOUNG_MODULE_Y][ii] * t.prop[YOUNG_MODULE][ii]; // ������ ���� ��.
				Ez = t.prop[MULT_YOUNG_MODULE_Z][ii] * t.prop[YOUNG_MODULE][ii]; // ������ ���� ��.

				nuyz = t.prop[MULT_POISSON_RATIO_YZ][ii] * t.prop[POISSON_RATIO][ii];
				nuxz = t.prop[MULT_POISSON_RATIO_XZ][ii] * t.prop[POISSON_RATIO][ii];
				nuxy = t.prop[MULT_POISSON_RATIO_XY][ii] * t.prop[POISSON_RATIO][ii];
				nuzy = t.prop[MULT_POISSON_RATIO_ZY][ii] * t.prop[POISSON_RATIO][ii];
				nuzx = t.prop[MULT_POISSON_RATIO_ZX][ii] * t.prop[POISSON_RATIO][ii];
				nuyx = t.prop[MULT_POISSON_RATIO_YX][ii] * t.prop[POISSON_RATIO][ii];
			}

			
			/*doublereal C_1[6][6] =
			{
				{1.0 / Ex, -nu / Ex, -nu / Ex, 0.0, 0.0, 0.0},
				{-nu / Ey, 1.0 / Ey, -nu / Ey, 0.0, 0.0, 0.0},
				{-nu / Ez, -nu / Ez, 1.0 / Ez, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0 / Gxy, 0.0, 0.0},
				{0.0, 0.0, 0.0, 0.0, 1.0 / Gyz, 0.0},
				{0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / Gxz}
			};*/

			
			doublereal C_1[6][6] =
			{
				{1.0 / Ex, -nuxy / Ex, -nuxz / Ex, 0.0, 0.0, 0.0},
				{-nuyx / Ey, 1.0 / Ey, -nuyz / Ey, 0.0, 0.0, 0.0},
				{-nuzx / Ez, -nuzy / Ez, 1.0 / Ez, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0 / Gxy, 0.0, 0.0},
				{0.0, 0.0, 0.0, 0.0, 1.0 / Gyz, 0.0},
				{0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / Gxz}
			};
			

			doublereal** Dirichlet = new doublereal * [6];
			for (int i_2 = 0; i_2 < 6; i_2++) {
				Dirichlet[i_2] = new doublereal[6];
				for (int i_3 = 0; i_3 < 6; i_3++) {
					Dirichlet[i_2][i_3] = C_1[i_2][i_3];
				}
			}

			inverse_matrix_simple(Dirichlet, 6, false);

			doublereal Strain_vec[6] = {
				t.total_deformation[STRAIN_X][i_1],
				t.total_deformation[STRAIN_Y][i_1],
				t.total_deformation[STRAIN_Z][i_1],
				t.total_deformation[STRAIN_XY][i_1],
				t.total_deformation[STRAIN_YZ][i_1],
				t.total_deformation[STRAIN_ZX][i_1]
			};

			doublereal Stress_vec[6] = {
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0
			};

			for (integer ir1 = 0; ir1 < 6; ir1++) {
				for (integer ir2 = 0; ir2 < 6; ir2++) {
					Stress_vec[ir1] += Dirichlet[ir1][ir2] * Strain_vec[ir2];
				}
			}

			t.total_deformation[STRESS_X][i_1] = Stress_vec[0];
			t.total_deformation[STRESS_Y][i_1] = Stress_vec[1];
			t.total_deformation[STRESS_Z][i_1] = Stress_vec[2];
			t.total_deformation[STRESS_XY][i_1] = Stress_vec[3];
			t.total_deformation[STRESS_YZ][i_1] = Stress_vec[4];
			t.total_deformation[STRESS_ZX][i_1] = Stress_vec[5];

			for (int i_2 = 0; i_2 < 6; i_2++) {
				delete[] Dirichlet[i_2];
			}
			delete[] Dirichlet;

			/*
			t.total_deformation[STRESS_X][i_1] = 0.0;
			for (integer i_11 = 0; i_11 < 6; i_11++) {
				if (((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
					// �� ������ ������ �� �������� epsilon_t ��� �� ���� ������������ � ������ �������������. 
					// � �� ������ ����� ����� ������ epsilon_t ��������� ??? 16.03.2021
					if ((i_11 == 0)) {
						t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] 
							/*- beta_t_solid_x * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
					if ((i_11 == 1)) {
						t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1]
							/*-	beta_t_solid_y * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
					if ((i_11 == 2)) {
						t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] 
							/*-	beta_t_solid_z * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
				}
				else {
					t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * t.total_deformation[i_11 + STRAIN_X][i_1];
				}
			}
			t.total_deformation[STRESS_Y][i_1] = 0.0;
			for (integer i_11 = 0; i_11 < 6; i_11++) {
				if (((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
					// �� ������ ������ �� �������� epsilon_t ��� �� ���� ������������ � ������ �������������. 
					// � �� ������ ����� ����� ������ epsilon_t ��������� ??? 16.03.2021
					if ((i_11 == 0)) {
						t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1]
							/*- beta_t_solid_x * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
					if ((i_11 == 1)) {
						t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1]
							/*- beta_t_solid_y * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
					if ((i_11 == 2)) {
						t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] 
							/*- beta_t_solid_z * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
				}
				else {
					t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * t.total_deformation[i_11 + STRAIN_X][i_1];
				}
			}
			t.total_deformation[STRESS_Z][i_1] = 0.0;
			for (integer i_11 = 0; i_11 < 6; i_11++) {
				if (((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
					// �� ������ ������ �� �������� epsilon_t ��� �� ���� ������������ � ������ �������������. 
					// � �� ������ ����� ����� ������ epsilon_t ��������� ??? 16.03.2021
					if ((i_11 == 0)) {
						t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] 
							/*-	beta_t_solid_x * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
					if ((i_11 == 1)) {
						t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] 
							/*-	beta_t_solid_y * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
					if ((i_11 == 2)) {
						t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1]
							/*-	beta_t_solid_z * (t.potent[i_1] - t.operatingtemperature_copy)*//*);
					}
				}
				else {
					t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * t.total_deformation[i_11 + STRAIN_X][i_1];
				}
			}

			t.total_deformation[STRESS_XY][i_1] = Dirichlet[3][3] * t.total_deformation[STRAIN_XY][i_1];
			t.total_deformation[STRESS_YZ][i_1] = Dirichlet[4][4] * t.total_deformation[STRAIN_YZ][i_1];
			t.total_deformation[STRESS_ZX][i_1] = Dirichlet[5][5] * t.total_deformation[STRAIN_ZX][i_1];
			*/
		}
		//for (integer i_11 = 0; i_11 < 6; i_11++) {
			//delete[] Dirichlet[i_11];
		//}
		//delete[] Dirichlet;

		/*
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxbound; i_1++) {
			t.total_deformation[STRESS_X][i_1] = t.total_deformation[STRESS_X][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_Y][i_1] = t.total_deformation[STRESS_Y][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_Z][i_1] = t.total_deformation[STRESS_Z][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_XY][i_1] = t.total_deformation[STRESS_XY][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_YZ][i_1] = t.total_deformation[STRESS_YZ][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_ZX][i_1] = t.total_deformation[STRESS_ZX][t.border_neighbor[i_1].iI];
		}
		*/
		// epsilon (STRAIN) von Mizes
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			//t.total_deformation[STRAIN_VON_MIZES][i_1] = sqrt(0.5*((t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Y][i_1])*
				//(t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Y][i_1]) + (t.total_deformation[STRAIN_Y][i_1] - t.total_deformation[STRAIN_Z][i_1]) *
				//(t.total_deformation[STRAIN_Y][i_1] - t.total_deformation[STRAIN_Z][i_1]) + (t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Z][i_1]) *
				//(t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Z][i_1])));

			t.total_deformation[STRAIN_VON_MIZES][i_1] = sqrt(0.5 * ((t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Y][i_1]) *
				(t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Y][i_1]) + (t.total_deformation[STRAIN_Y][i_1] - t.total_deformation[STRAIN_Z][i_1]) *
				(t.total_deformation[STRAIN_Y][i_1] - t.total_deformation[STRAIN_Z][i_1]) + (t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Z][i_1]) *
				(t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Z][i_1])+ 6.0 * (t.total_deformation[STRAIN_XY][i_1] * t.total_deformation[STRAIN_XY][i_1] +
					t.total_deformation[STRAIN_YZ][i_1] * t.total_deformation[STRAIN_YZ][i_1] + t.total_deformation[STRAIN_ZX][i_1] * t.total_deformation[STRAIN_ZX][i_1])));
					

			t.total_deformation[LOG10_STRAIN_VON_MIZES][i_1] = log10(t.total_deformation[STRAIN_VON_MIZES][i_1]);


			// STRESS

			
			//t.total_deformation[STRESS_VON_MIZES][i_1] = sqrt(0.5*((t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Y][i_1])*
				//(t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Y][i_1]) + (t.total_deformation[STRESS_Y][i_1] - t.total_deformation[STRESS_Z][i_1]) *
				//(t.total_deformation[STRESS_Y][i_1] - t.total_deformation[STRESS_Z][i_1]) + (t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Z][i_1]) *
				//(t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Z][i_1])));
				

			// https://ru.wikipedia.org/wiki/�����,_������_�����_���
			// https://ru.qaz.wiki/wiki/Von_Mises_yield_criterion
			t.total_deformation[STRESS_VON_MIZES][i_1] = sqrt(0.5 * ((t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Y][i_1]) *
				(t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Y][i_1]) + (t.total_deformation[STRESS_Y][i_1] - t.total_deformation[STRESS_Z][i_1]) *
				(t.total_deformation[STRESS_Y][i_1] - t.total_deformation[STRESS_Z][i_1]) + (t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Z][i_1]) *
				(t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Z][i_1]) + 6.0*(t.total_deformation[STRESS_XY][i_1] * t.total_deformation[STRESS_XY][i_1] +
				    t.total_deformation[STRESS_YZ][i_1] * t.total_deformation[STRESS_YZ][i_1]+ t.total_deformation[STRESS_ZX][i_1] * t.total_deformation[STRESS_ZX][i_1])));
					
			t.total_deformation[LOG10_STRESS_VON_MIZES][i_1] = log10(t.total_deformation[STRESS_VON_MIZES][i_1]);

		}

	
	}

	doublereal** total_deformation_shadow = nullptr;
	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
	{
		if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
			total_deformation_shadow = new doublereal*[20];
			for (integer j_1 = 0; j_1 < 20; j_1++) {
				total_deformation_shadow[j_1] = new doublereal[t.maxelm + t.maxbound];
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
					total_deformation_shadow[j_1][i_1] = t.total_deformation[j_1][i_1];
				}
			}
		}
	}


	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	int err = 0;
	switch (ikey) {
	case 0:  fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1:  fp = fopen64("ALICEFLOW0_27_temp.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
	default: fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
}
	if (fp==NULL) err = 1;
#else
	errno_t err;
	switch (ikey) {
	case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
	default: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	}
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2 in my_export_tecplot3.c\n");
		system("pause");

	}
	else {

		int c; // �������� ������
		integer ivarexport = 1; // �� ��������� ������ ���� ����������:
		integer i = 0; // ������� �����

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


			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count>0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // ������� ��� ����������� ������ �������
				}
			}

			if (ivarexport == 1) {

				integer ncell_shadow = 0;
#pragma omp parallel for reduction(+ : ncell_shadow)
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
					// ������������� ����������������� ��� ������� whot_is_block �����������
					// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
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
				
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// ������ ��� ����������
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
				if (1 && PHYSICAL_MODEL_SWITCH::MESHER_ONLY == steady_or_unsteady_global_determinant) {
					// ����� ������ ��������� ����������
					// ����������� ��������� network T �������.
					fprintf(fp, "VARIABLES = x, y, z, quolity, log(quolity), Temp \n");  // ���������� ������ ��������� �����
					//fprintf(fp, "VARIABLES = x, y, z, F(microW/cm!2), log(F(microW/cm!2)) \n");
				}
				else {
					// ���������� ������ ���� ����������
					if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
					{
						fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation, epsilon_x, epsilon_y, epsilon_z, gamma_xy, gamma_yz, gamma_zx, epsilon_von_Mizes, log10_epsilon_von_Mizes, stress_x, stress_y, stress_z, stress_xy, stress_yz, stress_zx, stress_von_Mizes, log10_stress_von_Mizes\n");  // ���������� ������ ���� ����������
					}
					else {
						if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
							fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux\n");  // ���������� ������ ���� ����������
						}
						else {
							fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");
						}
					}
				}

#if doubleintprecision == 1
				 // ������ ���������� � �����
				if (bTETRAHEDRON) {
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm + t.maxbound, 6*ncell_shadow);
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm, 6*ncell_shadow );
					}
				}
				else {
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow );
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow );
					}
				}
#else
					
				if (bTETRAHEDRON) {
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm + t.maxbound, 6 * ncell_shadow );
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm, 6 * ncell_shadow );
					}
				}
				else {
					// ������ ���������� � �����
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow );
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow );
					}
				}
#endif



				if (bvery_big_memory) {
					// extended printeger �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (1 && PHYSICAL_MODEL_SWITCH::MESHER_ONLY == steady_or_unsteady_global_determinant) {
						/*
						// ������ ������������.
						printf("calculation osveshennosti\n");
						system("PAUSE");
						doublereal* distance = new doublereal[t.maxelm + t.maxbound];
						for (i = 0; i < t.maxelm; i++) {
							distance[i] = 1.0e30;
							if (b[t.whot_is_block[i]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								distance[i] = 0.1e-3;// ���� ������� ��
							}
						}
						for (integer ibid = 0; ibid < lb; ibid++) {
							if (b[ibid].arr_Sc[0]>1.0e-30) {
								doublereal x0 = 0.5 * (b[ibid].g.xS + b[ibid].g.xE);
								doublereal y0 = 0.5 * (b[ibid].g.yS + b[ibid].g.yE);
								doublereal z0 = 0.5 * (b[ibid].g.zS + b[ibid].g.zE);
								for (i = 0; i < t.maxelm; i++) {
									TOCHKA p1;
									center_cord3D(i, t.nvtx, t.pa, p1, 100);
									doublereal dist0 = sqrt((p1.x-x0)*(p1.x-x0)+(p1.y-y0)*(p1.y-y0)+(p1.z-z0)*(p1.z-z0));
									if (dist0 < distance[i]) {
										distance[i] = dist0;
									}
								}
							}
						}

						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
							//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", 1.0/(M_PI* distance[i]* distance[i]));
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						*/

						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
							//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);

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
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
							//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
							volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
							doublereal dmax = dx;
							doublereal dmin = dx;
							if (dy > dmax) dmax = dy;
							if (dz > dmax) dmax = dz;
							if (dy < dmin) dmin = dy;
							if (dz < dmin) dmin = dz;

							fprintf(fp, "%+.6f ", log10(dmax / dmin));
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						// ����������� ��������� network T �������.
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.potent[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
				}
				else {
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
			}
			else if (ivarexport == 3) {
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ������ �� ������
					ncell_shadow = 0;

#pragma omp parallel for private(i) reduction(+: ncell_shadow)
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
					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

								integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
								integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
								integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
								integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


								integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
								integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
								integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
								integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



								integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
								integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
								integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
								integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];

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

								// ������������ ������ ������� ����.
								// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
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
										if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
										{
											for (integer j_4 = 0; j_4 < 20; j_4++) {
												total_deformation_shadow[j_4][inode5] = t.total_deformation[j_4][inode1];
												total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode2];
												total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode3];
												total_deformation_shadow[j_4][inode8] = t.total_deformation[j_4][inode4];
											}
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
										if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
										{
											for (integer j_4 = 0; j_4 < 20; j_4++) {
												total_deformation_shadow[j_4][inode2] = t.total_deformation[j_4][inode1];
												total_deformation_shadow[j_4][inode3] = t.total_deformation[j_4][inode4];
												total_deformation_shadow[j_4][inode6] = t.total_deformation[j_4][inode5];
												total_deformation_shadow[j_4][inode7] = t.total_deformation[j_4][inode8];
											}
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
										if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
											(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
										{
											for (integer j_4 = 0; j_4 < 20; j_4++) {
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

					}
					else {


						ncell_shadow = 0;

#pragma omp parallel for private(i) reduction(+: ncell_shadow)
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
					// � ��� ������ ���� ������� ����, ������� �� ���������� ������  ��������.
					ncell_shadow = ncell;
					ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;
				}
				if (1 && PHYSICAL_MODEL_SWITCH::MESHER_ONLY != steady_or_unsteady_global_determinant) {
					// ������ ����� ������� ������� � ������������� � �������������:

					if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
					{
						//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
						if ((f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
							fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, Re_y, k, omega, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation, epsilon_x, epsilon_y, epsilon_z, gamma_xy, gamma_yz, gamma_zx, epsilon_von_Mizes, log10_epsilon_von_Mizes, stress_x, stress_y, stress_z, stress_xy, stress_yz, stress_zx, stress_von_Mizes, log10_stress_von_Mizes\n");
						}
						else {
							fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, Re_y, k, eps, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation, epsilon_x, epsilon_y, epsilon_z, gamma_xy, gamma_yz, gamma_zx, epsilon_von_Mizes, log10_epsilon_von_Mizes, stress_x, stress_y, stress_z, stress_xy, stress_yz, stress_zx, stress_von_Mizes, log10_stress_von_Mizes\n");
						}
					}
					else {
						if ((f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) || (f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
							fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, Re_y, k, omega, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux\n");
						}
						else {
							fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, Re_y, k, eps, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux\n");
						}
					}
				}
				else {
					// 19,03,2019
					//fprintf(fp, "\nVARIABLES = x, y, z, quolity, log(quolity)\n"); // ������ ����� ������
					fprintf(fp, "\nVARIABLES = x, y, z, F(microW/cm!2), log(F(microW/cm!2))\n");
				}

#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bTETRAHEDRON) {
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm + t.maxbound, 6 * ncell_shadow);
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm, 6 * ncell_shadow );
					}
				}
				else {
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow );
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow );
					}
				}
#else
				if (bTETRAHEDRON) {
					// ������ ���������� � �����
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm + t.maxbound, 6*ncell_shadow );
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TETRAHEDRON, F=FEBLOCK\n\n", maxelm, 6*ncell_shadow );
					}
				}
				else {
					// ������ ���������� � �����
					if (bextendedprint) {
						fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow );
					}
					else {
						fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow );
					}
			}
#endif




				if (bvery_big_memory) {
					if (lite_export) {
						// ������ x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						if (1 && (PHYSICAL_MODEL_SWITCH::MESHER_ONLY == steady_or_unsteady_global_determinant)) {
							/*
													// ������ ������������.
							printf("calculation osveshennosti\n");
							system("PAUSE");
							doublereal* distance = new doublereal[t.maxelm + t.maxbound];
							for (i = 0; i < t.maxelm; i++) {
								distance[i] = 1.0e30;
								if (b[t.whot_is_block[i]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
									distance[i] = -0.1e-3;// ���� ������� ��
								}
							}
							for (integer ibid = 0; ibid < lb; ibid++) {
								if (b[ibid].arr_Sc[0] > 1.0e-30) {
									printf("power %e\n", b[ibid].arr_Sc[0]);
									doublereal x0 = 0.5 * (b[ibid].g.xS + b[ibid].g.xE);
									doublereal y0 = 0.5 * (b[ibid].g.yS + b[ibid].g.yE);
									doublereal z0 = 0.5 * (b[ibid].g.zS + b[ibid].g.zE);
									for (i = 0; i < t.maxelm; i++) {
										TOCHKA p1;
										center_cord3D(i, t.nvtx, t.pa, p1, 100);
										doublereal dist0 = sqrt((p1.x - x0) * (p1.x - x0) + (p1.y - y0) * (p1.y - y0) + (p1.z - z0) * (p1.z - z0));
										if (dist0 < distance[i]) {
											distance[i] = dist0;
										}
									}
								}
							}
							for (i = 0; i < t.maxelm; i++) {
								if (distance[i] < 0.0) {
									distance[i] = 1.0e10;
								}
							}
							//for (i = 0; i < t.maxelm; i++) {

							//}
							//for (i = 0; i < t.maxbound; i++) {

							//}

							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
								//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.6f ", 1.0 / (M_PI * distance[i] * distance[i]));
								if (i % 10 == 0) fprintf(fp, "\n");
							}

							delete[] distance;
							*/
							// ������ ������������.
							printf("calculation illumination\n");
							system("PAUSE");
							doublereal* myF = new doublereal[t.maxelm + t.maxbound];

							calculate_light_flux(myF, t, b, lb);

							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
								//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.6f ", myF[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}



							/*
														// quolity
														for (i = 0; i < t.database.maxelm; i++) {
															doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
															//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
															volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
															doublereal dmax = dx;
															doublereal dmin = dx;
															if (dy > dmax) dmax = dy;
															if (dz > dmax) dmax = dz;
															if (dy < dmin) dmin = dy;
															if (dz < dmin) dmin = dz;

															fprintf(fp, "%+.6f ", dmax/dmin);
															if (i % 10 == 0) fprintf(fp, "\n");
														}
														*/
														// log10(quolity)
							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
								//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								//fprintf(fp, "%+.6f ", log10(dmax / dmin));
								fprintf(fp, "%+.6f ", fmax(0.0, log10(myF[i])));
								if (i % 10 == 0) fprintf(fp, "\n");
							}

							delete[] myF;

							// ����������� ��������� network T �������.
							///for (i = 0; i < t.database.maxelm; i++) {
								//fprintf(fp, "%+.16f ", t.potent[i]);
								//if (i % 10 == 0) fprintf(fp, "\n");
							//}

							fprintf(fp, "\n");

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

									if (bTETRAHEDRON) {

										printf("incomming");
										system("pause");

										print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
											b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
											t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
											t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
									}
									else {
										if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
											// ������������ ������� ���� � ������ ��� � ������.
											//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
											//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);

											fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);


#else
											// ������������ ������� ���� � ������ ��� � ������.
											//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
											//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
											fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


										}
									}
								}
							}
						}
					}
					else {
						// ������ x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						// ������ z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}

						if (1 && (PHYSICAL_MODEL_SWITCH::MESHER_ONLY == steady_or_unsteady_global_determinant)) {
							// quolity
							for (i = 0; i < t.database.maxelm; i++) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
								//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
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
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
								//volume3D(i, t.nvtx, t.pa, dx, dy, dz);
								volume3D_q(i, t.nvtx, t.pa, dx, dy, dz);
								doublereal dmax = dx;
								doublereal dmin = dx;
								if (dy > dmax) dmax = dy;
								if (dz > dmax) dmax = dz;
								if (dy < dmin) dmin = dy;
								if (dz < dmin) dmin = dz;

								fprintf(fp, "%+.16f ", log10(dmax / dmin));
								if (i % 10 == 0) fprintf(fp, "\n");
							}

							// ����������� ��������� network T �������.
							for (i = 0; i < t.database.maxelm; i++) {								
								fprintf(fp, "%+.16f ", t.potent[i]);
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
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// ������ ������ �����

		// ������ ���� ���������� ������������ ������.
		if (1 && PHYSICAL_MODEL_SWITCH::MESHER_ONLY != steady_or_unsteady_global_determinant) {

		// ������ �����������
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%+.3f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.3f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%+.16f ", temp_shadow[i]);
				}
				else {
					fprintf(fp, "%+.16f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx*svx + svy*svy + svz*svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
					doublereal svy = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
					doublereal svz = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // ������������ ������������ ��������
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);

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
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm] / f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if ((f[0].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[0].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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
					if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
						doublereal speed_or_sqrt_k = 0.0;
						if (f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][t.ptr[0][i]]));
						}
						if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
							(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY][t.ptr[0][i]]));
						}
						//speed_or_sqrt_k = sqrt(f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] + f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] + f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
						integer id_loc = t.ptr[0][i];
						doublereal Re_y = 0.0;
						if (id_loc < f[t.ptr[1][i]].maxelm) {
							Re_y = speed_or_sqrt_k *
								f[t.ptr[1][i]].rdistWall[t.ptr[0][i]] * f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]];
						}
						else {
							Re_y = speed_or_sqrt_k *
								f[t.ptr[1][i]].rdistWall[t.ptr[0][i]] * f[t.ptr[1][i]].prop_b[RHO][t.ptr[0][i] - f[t.ptr[1][i]].maxelm] / f[t.ptr[1][i]].prop_b[MU_DYNAMIC_VISCOSITY][t.ptr[0][i] - f[t.ptr[1][i]].maxelm];

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
					if ((f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
						doublereal speed_or_sqrt_k = 0.0;
						if (f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[idfluid].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i + maxelm]));
						}
						if ((f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
							(f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
							speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f[idfluid].potent[TURBULENT_KINETIK_ENERGY][i + maxelm]));
						}
						//speed_or_sqrt_k = sqrt(f[idfluid].potent[VELOCITY_X_COMPONENT][i+ maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i+ maxelm] + f[idfluid].potent[VELOCITY_Y_COMPONENT][i+ maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i+ maxelm] + f[idfluid].potent[VELOCITY_Z_COMPONENT][i+ maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i+ maxelm]);
						doublereal Re_y = speed_or_sqrt_k *
							f[idfluid].rdistWall[i + maxelm] * f[idfluid].prop_b[RHO][i] / f[idfluid].prop_b[MU_DYNAMIC_VISCOSITY][i];

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
					if (f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][t.ptr[0][i]]);
					}
					else if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if (f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i + maxelm]); // TURBULENT_KINETIK_ENERGY_STD_K_EPS
					}
					else if ((f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if (f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][t.ptr[0][i]]);
					}
					else if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if (f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i + maxelm]); // TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS
					}
					else if ((f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
						(f[idfluid].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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

		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {

		doublereal *Tx = nullptr;
		doublereal *Ty = nullptr;
		doublereal *Tz = nullptr;
		Tx = new doublereal[t.maxelm + t.maxbound];
		Ty = new doublereal[t.maxelm + t.maxbound];
		Tz = new doublereal[t.maxelm + t.maxbound];

		// ������������� ����.
#pragma omp parallel for
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// ���������� ����������.
#pragma omp parallel for private(i)
		for (i = 0; i<t.maxelm; i++) {
			// ������ ���������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
		}

#pragma omp parallel for private(i)
		for (i = 0; i<t.maxelm; i++) {
			// ������ ��������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
		}


		// ���������� � ����.

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


		// ������������ ����������� ������.
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

		fprintf(fp, "\n");

		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{

			doublereal dmax_d = -1.0e30;
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					// total deformation
					doublereal t1 =
						sqrt(total_deformation_shadow[1][i] * total_deformation_shadow[1][i] +
							total_deformation_shadow[2][i] * total_deformation_shadow[2][i] +
							total_deformation_shadow[3][i] * total_deformation_shadow[3][i]);

					if (fabs(t1) > 1.0e-20) {
						if (total_deformation_shadow[0][i] / t1 > dmax_d) {
							dmax_d = total_deformation_shadow[0][i] / t1;
						}
					}
					total_deformation_shadow[0][i] = t1;
				}
				else {
					// total deformation
					doublereal t1 =
						sqrt(t.total_deformation[1][i] * t.total_deformation[1][i] +
							t.total_deformation[2][i] * t.total_deformation[2][i] +
							t.total_deformation[3][i] * t.total_deformation[3][i]);

					if (fabs(t1) > 1.0e-20) {
						if (t.total_deformation[0][i] / t1 > dmax_d) {
							dmax_d = t.total_deformation[0][i] / t1;
						}
					}
					t.total_deformation[0][i] = t1;
				}
			}
			printf("otklonenie =%e\n", dmax_d);

		for (integer j_6 = 0; j_6 < 20; j_6++) {

			// ������ ������ ����������
			if (lite_export) {
				for (i = 0; i < maxelm; i++) {


					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
						fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
					}
					else {
						fprintf(fp, "%e ", t.total_deformation[j_6][i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
		} // if (1 && MESHER_ONLY != steady_or_unsteady_global_determinant) 
		}


		if (bvery_big_memory) {
			// ������ ���������� � ���������� �����
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

					if (bTETRAHEDRON) {

						printf("incomming");
						system("pause");

						print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
							b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
							t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
							t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
					}
					else {
						if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
							// ������������ ������� ���� � ������ ��� � ������.
							//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							
							
#else
							// ������������ ������� ���� � ������ ��� � ������.
							//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}
					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
					{

						// �� �������� � TETRAHEDRON

						integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
						inode1 = t.database.nvtxcell[0][i] - 1;
						inode2 = t.database.nvtxcell[1][i] - 1;
						inode3 = t.database.nvtxcell[2][i] - 1;
						inode4 = t.database.nvtxcell[3][i] - 1;
						inode5 = t.database.nvtxcell[4][i] - 1;
						inode6 = t.database.nvtxcell[5][i] - 1;
						inode7 = t.database.nvtxcell[6][i] - 1;
						inode8 = t.database.nvtxcell[7][i] - 1;

						integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
						integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
						integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
						integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


						integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
						integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
						integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
						integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



						integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
						integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
						integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
						integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];


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

						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if (bTETRAHEDRON) {

								print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
									b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
									t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
									t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							}
							else {

								if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

								}
							}
						}
						else if (((inode5B >= 0) && (inode5B < t.maxelm) && (inode6B >= 0) && (inode6B < t.maxelm) && (inode7B >= 0) && (inode7B < t.maxelm) && (inode8B >= 0) && (inode8B < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (!((t.ptr[1][inode5B] == -1) && (t.ptr[1][inode6B] == -1) && (t.ptr[1][inode7B] == -1) && (t.ptr[1][inode8B] == -1))) && (!((t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

								if (bTETRAHEDRON) {

									print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
										b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
										t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
										t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								}
								else {
#if doubleintprecision == 1
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif
								}
							}
						}
						else if (((inode2W >= 0) && (inode2W < t.maxelm) && (inode3W >= 0) && (inode3W < t.maxelm) && (inode6W >= 0) && (inode6W < t.maxelm) && (inode7W >= 0) && (inode7W < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode8] == -1) && (!((t.ptr[1][inode2W] == -1) && (t.ptr[1][inode3W] == -1) && (t.ptr[1][inode6W] == -1) && (t.ptr[1][inode7W] == -1))) && (!((t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1)))))
						{
							if ((b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible)) {

								if (bTETRAHEDRON) {

									print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
										b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
										t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
										t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								}
								else {

#if doubleintprecision == 1
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif
								}
							}
						}
						else if (((inode3S >= 0) && (inode3S < t.maxelm) && (inode4S >= 0) && (inode4S < t.maxelm) && (inode7S >= 0) && (inode7S < t.maxelm) && (inode8S >= 0) && (inode8S < t.maxelm) && (t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (!((t.ptr[1][inode3S] == -1) && (t.ptr[1][inode4S] == -1) && (t.ptr[1][inode7S] == -1) && (t.ptr[1][inode8S] == -1))) && (!((t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))))
						{
							if ((b[ib4].bvisible) && (b[ib3].bvisible) && (b[ib8].bvisible) && (b[ib7].bvisible)) {

								if (bTETRAHEDRON) {

									print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
										b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
										t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
										t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								}
								else {
#if doubleintprecision == 1
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
									fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif
								}
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

						if (bTETRAHEDRON) {

							print_tetrahedron(fp, b[ib1].bvisible, b[ib2].bvisible, b[ib3].bvisible, b[ib4].bvisible,
								b[ib5].bvisible, b[ib6].bvisible, b[ib7].bvisible, b[ib8].bvisible,
								t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i],
								t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						}
						else {

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								// ������������ ������� ���� � ��������.
								// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
								// ������������ ������� ���� � ��������.
								// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
	}

	if (temp_shadow != nullptr) {
		delete[] temp_shadow;
		temp_shadow = nullptr;
	}
	//total_deformation
	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
	{
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
	}

	

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);



} // exporttecplotxy360T_3D_part2

void export_tecplot_temperature_ass(int** &nvtx, bool* &bcheck_visible, TOCHKA* &pa, doublereal* &potent,
	doublereal* &lam_for_export, doublereal* &Txgl, doublereal* &Tygl, doublereal* &Tzgl, doublereal* &HeatFluxMaggl,
	 integer maxelm, integer ncell) {
	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp=fopen64("ALICEFLOW0_08_temp.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_08_temp.PLT", "wb");
#endif
	

	if ((err) != 0) {
		printf("Create File temp Error in function export_tecplot_temperature_ass in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"ALICEFLOW0_08\"\n");

		// ������ ��� ����������
		fprintf(fp, "VARIABLES = x, y, z, T, lam, Tx, Ty, Tz, magHeatFlux, log10MagHeatFlux \n");  // ���������� ������ ���� ����������

												   // ������ ���������� � �����
		integer ncell_actual = 0;
		for (integer i = 0; i < ncell; i++) {
			if (bcheck_visible[i]) {
				ncell_actual++;
			}
		}
		fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_actual);


		// ������ x
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", pa[i].x);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ y
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", pa[i].y);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ z
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", pa[i].z);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ �����������
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", potent[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ ����������������
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", lam_for_export[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		
		// ������ -lambda*gradTx
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", Txgl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ -lambda*gradTy
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", Tygl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ -lambda*gradTz
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", Tzgl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ heatfluxMag
		for (integer i = 0; i < maxelm; i++) {
			fprintf(fp, "%+.6f ", HeatFluxMaggl[i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ log10_heatfluxMag
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
				// ��� �� ���� �� ���� ��������� ������� � nvtx. ����� �� ������������.
				if (b_on_adaptive_local_refinement_mesh) {
					fprintf(fp, "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
				}
				else {
					fprintf(fp, "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[3][i], nvtx[2][i], nvtx[4][i], nvtx[5][i], nvtx[7][i], nvtx[6][i]);
				}
			}
		}

		fclose(fp);
		printf("file succsefull writing\n");

	}
}


void exporttecplot_assembles_mesh(TEMPER &t, integer lu, UNION* &my_union) {
	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");
#endif
	

	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplot_assembles_mesh in my_export_tecplot3.c\n");
		system("pause");

	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

		// ������ ��� ����������
		fprintf(fp, "VARIABLES = x, y, z \n");  // ���������� ������ ���� ����������

		

		integer maxelm_global = t.maxnod;
		integer ncell_shadow_gl = t.maxelm;
		for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
			maxelm_global += my_union[iu_74].t.maxnod;			
			ncell_shadow_gl += my_union[iu_74].t.maxelm;
		}


		// ������ ���������� � �����
		
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

			// ������ x
			for (integer i = 0; i < t.maxnod; i++) {
				fprintf(fp, "%+.6f ", t.pa[i].x);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// ������ x
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					fprintf(fp, "%+.6f ", my_union[iunion_scan].t.pa[i].x);
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}
			

			fprintf(fp, "\n");

			// ������ y
			for (integer i = 0; i < t.maxnod; i++) {
				fprintf(fp, "%+.6f ", t.pa[i].y);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// ������ x
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					fprintf(fp, "%+.6f ", my_union[iunion_scan].t.pa[i].y);
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}
			

			fprintf(fp, "\n");

			// ������ z
			for (integer i = 0; i < t.maxnod; i++) {
				fprintf(fp, "%+.6f ", t.pa[i].z);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// ������ z
				for (integer i = 0; i < my_union[iunion_scan].t.maxnod; i++) {
					fprintf(fp, "%+.6f ", my_union[iunion_scan].t.pa[i].z);
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}
			
			fprintf(fp, "\n");

			int iadd = 0;

			for (integer i = 0; i < t.maxelm; i++) {
				fprintf(fp, "%d %d %d %d %d %d %d %d \n", t.nvtx[0][i], t.nvtx[1][i], t.nvtx[3][i], t.nvtx[2][i], t.nvtx[4][i], t.nvtx[5][i], t.nvtx[7][i], t.nvtx[6][i]);
			}
			
			iadd += t.maxnod;

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				for (integer i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					fprintf(fp, "%d %d %d %d %d %d %d %d \n", iadd+my_union[iunion_scan].t.nvtx[0][i], iadd + my_union[iunion_scan].t.nvtx[1][i], iadd + my_union[iunion_scan].t.nvtx[3][i], iadd + my_union[iunion_scan].t.nvtx[2][i], iadd + my_union[iunion_scan].t.nvtx[4][i], iadd + my_union[iunion_scan].t.nvtx[5][i], iadd + my_union[iunion_scan].t.nvtx[7][i], iadd + my_union[iunion_scan].t.nvtx[6][i]);
				}
				iadd += my_union[iunion_scan].t.maxnod;
			}

			if (pa_gl != nullptr) {
				delete[] pa_gl;
			}

			fclose(fp);

	}

}

  // 10 ������ 2016 . �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
  // � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
  // �������� ����������� �����
  // ������� ���������� ������� � ��������� tecplot360
  // ����� 2.
void exporttecplotxy360T_3D_part2_assembles(integer maxelm, integer ncell, 
	FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, 
	bool bextendedprint, integer ikey, integer lu, UNION* &my_union,
	BLOCK* &b, integer &lb)
{
	const bool lite_export = true;
	// 16 ������ ����� ������� ��������� ������ �������,
	// ������ ���������� ����� ������.
	switch (ionly_solid_visible) {
	case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE:
		std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
		break;
	case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE:
		std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
		break;
	}
	

	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // ����� 1 ��� 3
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	doublereal** temp_shadow = nullptr;
	temp_shadow = new doublereal*[lu+1];
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
		temp_shadow[0] = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[0][i_1] = t.potent[i_1];
		}
	}
	for (integer i = 0; i < lu; i++) {
		if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
			temp_shadow[i+1] = new doublereal[my_union[i].t.maxelm + my_union[i].t.maxbound];
			for (integer i_1 = 0; i_1 < my_union[i].t.maxelm + my_union[i].t.maxbound; i_1++) {
				temp_shadow[i+1][i_1] = my_union[i].t.potent[i_1];
			}
		}
	}

	doublereal*** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
		if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
			total_deformation_shadow[i+1] = new doublereal*[4];
			for (integer j_1 = 0; j_1 < 4; j_1++) {
				total_deformation_shadow[i+1][j_1] = new doublereal[my_union[i].t.maxelm + my_union[i].t.maxbound];
				for (integer i_1 = 0; i_1 < my_union[i].t.maxelm + my_union[i].t.maxbound; i_1++) {
					total_deformation_shadow[i+1][j_1][i_1] = my_union[i].t.total_deformation[j_1][i_1];
				}
			}
		}
	}

	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	int err = 0;
	switch (ikey) {
	  case 0: fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
	  case 1: fp = fopen64("ALICEFLOW0_27_temp.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
	  default: fp = fopen64("ALICEFLOW0_07_temp.PLT", "wb");  break;
    }
	if (fp==NULL) err = 1;
#else
	errno_t err;
	switch (ikey) {
	case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
	default: err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "wb");  break;
	}
#endif
	
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2_assembles in file my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // �������� ������
		integer ivarexport = 1; // �� ��������� ������ ���� ����������:
		integer i = 0; // ������� �����

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


			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count>0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // ������� ��� ����������� ������ �������
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
					// ������������� ����������������� ���-������� whot_is_block �����������
					// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� ������.
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
						// ������������� ����������������� ���-������� whot_is_block �����������
						// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� ������.
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
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// ������ ��� ����������
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
				fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // ���������� ������ ���� ����������

				integer maxelm_global = maxelm;
				integer maxbound_global = t.maxbound;
				integer ncell_shadow_gl = ncell_shadow[0];
				for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
					maxelm_global += my_union[iu_74].t.maxelm;
					maxbound_global += my_union[iu_74].t.maxbound;
					ncell_shadow_gl += ncell_shadow[iu_74+1];
				}

#if doubleintprecision == 1
																																																							   // ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#else
																																																							   // ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#endif



				if (bvery_big_memory) {
					// extended print �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
						// ������ x
						for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
						// ������ y
						for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
					}
					// ������ z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
						// ������ z
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
				// ������ ���������
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
					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

							integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
							integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
							integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
							integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


							integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
							integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
							integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
							integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



							integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
							integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
							integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
							integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];

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

							// ������������ ������ ������� ����.
							// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
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

								integer inode2W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
								integer inode3W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
								integer inode6W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
								integer inode7W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


								integer inode5B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
								integer inode6B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
								integer inode7B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
								integer inode8B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



								integer inode3S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
								integer inode4S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
								integer inode7S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
								integer inode8S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode5];

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

								// ������������ ������ ������� ����.
								// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
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
							// � ��� ������ ���� ������� ����, ������� �� ���������� ������  ��������.
							ncell_shadow[0] = ncell;
							for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
								ncell_shadow[iunion_scan+1] += my_union[iunion_scan].t.ncell;
							}
							ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;
						}
					}
				}
				// ������ ����� ������� ������� � ������������� � �������������:
				
				//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
				fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
				

				integer maxelm_global = maxelm;
				integer maxbound_global = t.maxbound;
				integer ncell_shadow_gl = ncell_shadow[0];
				for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
					maxelm_global += my_union[iu_74].t.maxelm;
					maxbound_global += my_union[iu_74].t.maxbound;
					ncell_shadow_gl+= ncell_shadow[iu_74+1];
				}

#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global + maxbound_global, ncell_shadow_gl);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm_global, ncell_shadow_gl);
				}
#endif

				if (bvery_big_memory) {
					if (lite_export) {
						// ������ x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// ������ x
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", my_union[iunion_scan].t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// ������ y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// ������ y
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", my_union[iunion_scan].t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// ������ z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.6f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// ������ z
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", my_union[iunion_scan].t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
					else {
						// ������ x
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.x[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// ������ x
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// ������ y
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.y[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// ������ y
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						// ������ z
						for (i = 0; i < t.database.maxelm; i++) {
							fprintf(fp, "%+.16f ", t.database.z[i]);
							if (i % 10 == 0) fprintf(fp, "\n");
						}
						for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
							// ������ z
							for (i = 0; i < my_union[iunion_scan].t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", my_union[iunion_scan].t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
					}
				}
				else {
					// ��� ����� ���� ������ ������� �� ������������,
					// ��� �������� � ����������� ������, �.�. ������ � 
					// ���������� ����� ������� ���������.
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}

				delete[] ncell_shadow;
			}
			if (!bvery_big_memory) {
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// ������ ������ �����

		// ������ ���� ���������� ������������ ������.


		// ������ �����������
		if (lite_export) {
			for (i = 0; i < maxelm; i++) {
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%+.3f ", temp_shadow[0][i]);
				}
				else {
					fprintf(fp, "%+.3f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
						fprintf(fp, "%+.3f ", temp_shadow[iunion_scan + 1][i]);
					}
					else {
						fprintf(fp, "%+.3f ", my_union[iunion_scan].t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = my_union[iunion_scan].t.maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
				if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
					fprintf(fp, "%+.16f ", temp_shadow[0][i]);
				}
				else {
					fprintf(fp, "%+.16f ", t.potent[i]);
				}
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = maxelm; i < maxelm + t.maxbound; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
						fprintf(fp, "%+.16f ", temp_shadow[iunion_scan + 1][i]);
					}
					else {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = my_union[iunion_scan].t.maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
					doublereal svy = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
					doublereal svz = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Speed
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						doublereal svx = my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][my_union[iunion_scan].t.ptr[0][i]] * my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][my_union[iunion_scan].t.ptr[0][i]];
						doublereal svy = my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][my_union[iunion_scan].t.ptr[0][i]] * my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][my_union[iunion_scan].t.ptr[0][i]];
						doublereal svz = my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][my_union[iunion_scan].t.ptr[0][i]] * my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][my_union[iunion_scan].t.ptr[0][i]];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}


				if (bextendedprint) {
					// Speed
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						doublereal svx = my_union[iunion_scan].f[idfluid].potent[VELOCITY_X_COMPONENT][i + my_union[iunion_scan].t.maxelm] * my_union[iunion_scan].f[idfluid].potent[VELOCITY_X_COMPONENT][i + my_union[iunion_scan].t.maxelm];
						doublereal svy = my_union[iunion_scan].f[idfluid].potent[VELOCITY_Y_COMPONENT][i + my_union[iunion_scan].t.maxelm] * my_union[iunion_scan].f[idfluid].potent[VELOCITY_Y_COMPONENT][i + my_union[iunion_scan].t.maxelm];
						doublereal svz = my_union[iunion_scan].f[idfluid].potent[VELOCITY_Z_COMPONENT][i + my_union[iunion_scan].t.maxelm] * my_union[iunion_scan].f[idfluid].potent[VELOCITY_Z_COMPONENT][i + my_union[iunion_scan].t.maxelm];
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// VX
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VX
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[VELOCITY_X_COMPONENT][i + my_union[iunion_scan].t.maxelm]); // VX
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// VY
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VY
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[VELOCITY_Y_COMPONENT][i + my_union[iunion_scan].t.maxelm]); // VY
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// VZ
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][my_union[iunion_scan].t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VZ
					integer idfluid = 0;
					for (i = 0; i < my_union[iunion_scan].f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[VELOCITY_Z_COMPONENT][i + my_union[iunion_scan].t.maxelm]); // VZ
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Mu
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][my_union[iunion_scan].t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Mu
					for (i = 0; i < my_union[iunion_scan].f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");

			// Mut // ������������ ������������ ��������
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);

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
					fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm] / f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// Mut // ������������ ������������ ��������
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[MUT][my_union[iunion_scan].t.ptr[0][i]]);
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].potent[MUT][my_union[iunion_scan].t.ptr[0][i]] / my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][my_union[iunion_scan].t.ptr[0][i]]);

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
						fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].potent[MUT][i + my_union[iunion_scan].t.maxelm] / my_union[iunion_scan].f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

			}

			fprintf(fp, "\n");

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if ((f[0].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[0].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				// ��� ������� ������ ������������� ��������� ����������� �������.
				// ��� ������������� ���������� �� ������ Distance_Wall.
				for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {
					if (my_union[iunion_scan].t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", doublereal(i));
						if ((my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
							(my_union[iunion_scan].f[my_union[iunion_scan].t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
						if ((my_union[iunion_scan].f[0].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
							(my_union[iunion_scan].f[0].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
							(my_union[iunion_scan].f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES)||
							(my_union[iunion_scan].f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
							(my_union[iunion_scan].f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
							(my_union[iunion_scan].f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
							fprintf(fp, "%+.16f ", my_union[iunion_scan].f[idfluid].rdistWall[i + my_union[iunion_scan].t.maxelm]); // Distance_Wall
						}
						else fprintf(fp, "%+.16f ", 0.0);
						if ((i + my_union[iunion_scan].t.maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}

			fprintf(fp, "\n");



			// Curl // ������������ - ������ ������ ��������.
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
				// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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

		// ������������� ����.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[0][i] = 0.0;
			Ty[0][i] = 0.0;
			Tz[0][i] = 0.0;
		}

		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// ������������� ����.
			for (i = 0; i<my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
				Tx[iunion_scan + 1][i] = 0.0;
				Ty[iunion_scan + 1][i] = 0.0;
				Tz[iunion_scan + 1][i] = 0.0;
			}
		}

		// ���������� ����������.
		for (i = 0; i<t.maxelm; i++) {
			// ������ ���������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, Tx[0], Ty[0], Tz[0], t.ilevel_alice);
		}

		for (i = 0; i<t.maxelm; i++) {
			// ������ ��������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, Tx[0], Ty[0], Tz[0], t.ilevel_alice);
		}


		for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
			// ���������� ����������.
			for (i = 0; i<my_union[iunion_scan].t.maxelm; i++) {
				// ������ ���������� ����.
				green_gaussTemperature(i, my_union[iunion_scan].t.potent, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa,
					my_union[iunion_scan].t.neighbors_for_the_internal_node, my_union[iunion_scan].t.maxelm, false,
					my_union[iunion_scan].t.border_neighbor, Tx[iunion_scan + 1], Ty[iunion_scan + 1], Tz[iunion_scan + 1], 
					my_union[iunion_scan].t.ilevel_alice);
			}

			for (i = 0; i<my_union[iunion_scan].t.maxelm; i++) {
				// ������ ��������� ����.
				green_gaussTemperature(i, my_union[iunion_scan].t.potent, my_union[iunion_scan].t.nvtx, my_union[iunion_scan].t.pa,
					my_union[iunion_scan].t.neighbors_for_the_internal_node, my_union[iunion_scan].t.maxelm, true,
					my_union[iunion_scan].t.border_neighbor, Tx[iunion_scan + 1], Ty[iunion_scan + 1], Tz[iunion_scan + 1],
					my_union[iunion_scan].t.ilevel_alice);
			}
		}

		// ���������� � ����.

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
			//const doublereal eps_min = 2.0;
			//const doublereal d_stub = 0.0;
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
						i = i1;
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


		// ������������ ����������� ������.

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

			// ������ ������ ����������
			if (lite_export) {
				for (i = 0; i < maxelm; i++) {


					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
						fprintf(fp, "%e ", total_deformation_shadow[0][j_6][i]);
					}
					else {
						fprintf(fp, "%e ", t.total_deformation[j_6][i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
				// ������ ������ ����������
				if (lite_export) {
					for (i = 0; i < my_union[iunion_scan].t.maxelm; i++) {


						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
							if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
							fprintf(fp, "%e ", total_deformation_shadow[iunion_scan + 1][j_6][i]);
						}
						else {
							fprintf(fp, "%e ", my_union[iunion_scan].t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (bextendedprint) {
						for (i = maxelm; i < my_union[iunion_scan].t.maxelm + my_union[iunion_scan].t.maxbound; i++) {
							if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
			// ������ ���������� � ���������� �����
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
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

						integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
						integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
						integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
						integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


						integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
						integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
						integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
						integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



						integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
						integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
						integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
						integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];


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

						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
							// ������������ ������� ���� � ��������.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// ������������ ������� ���� � ��������.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}

			// TODO
			// ������������ ������� ��������� ���������� �����.
			int iadd = t.maxelm;
			for (integer iunion_scan = 0; iunion_scan < lu; iunion_scan++) {
				if (iunion_scan > 0) {
					iadd += my_union[iunion_scan-1].t.maxelm;
				}

				// ������ ���������� � ���������� �����
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
							// ������������ ������� ���� � ������ ��� � ������.
							//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i], iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
							// ������������ ������� ���� � ������ ��� � ������.
							//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i], iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
#endif


						}

					}
					else {
						//printf("fluid plot\n");
						//getchar();

						if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

							integer inode2W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
							integer inode3W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
							integer inode6W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
							integer inode7W = my_union[iunion_scan].t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


							integer inode5B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
							integer inode6B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
							integer inode7B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
							integer inode8B = my_union[iunion_scan].t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



							integer inode3S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
							integer inode4S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
							integer inode7S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
							integer inode8S = my_union[iunion_scan].t.neighbors_for_the_internal_node[S_SIDE][0][inode5];


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

							// ������������ ������ ������� ����.
							// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
							if (((my_union[iunion_scan].t.ptr[1][inode1] == -1) && (my_union[iunion_scan].t.ptr[1][inode2] == -1) && (my_union[iunion_scan].t.ptr[1][inode3] == -1) && (my_union[iunion_scan].t.ptr[1][inode4] == -1) &&
								(my_union[iunion_scan].t.ptr[1][inode5] == -1) && (my_union[iunion_scan].t.ptr[1][inode6] == -1) && (my_union[iunion_scan].t.ptr[1][inode7] == -1) && (my_union[iunion_scan].t.ptr[1][inode8] == -1)))
							{

								if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
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
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
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
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
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
									fprintf(fp, "%d %d %d %d ", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd + my_union[iunion_scan].t.database.nvtxcell[1][i], iadd + my_union[iunion_scan].t.database.nvtxcell[2][i], iadd + my_union[iunion_scan].t.database.nvtxcell[3][i]);
									fprintf(fp, "%d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[4][i], iadd + my_union[iunion_scan].t.database.nvtxcell[5][i], iadd + my_union[iunion_scan].t.database.nvtxcell[6][i], iadd + my_union[iunion_scan].t.database.nvtxcell[7][i]);
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
								// ������������ ������� ���� � ��������.
								// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
								fprintf(fp, "%d %d %d %d %d %d %d %d\n", iadd + my_union[iunion_scan].t.database.nvtxcell[0][i], iadd+my_union[iunion_scan].t.database.nvtxcell[1][i], iadd+my_union[iunion_scan].t.database.nvtxcell[2][i], iadd+my_union[iunion_scan].t.database.nvtxcell[3][i], iadd+my_union[iunion_scan].t.database.nvtxcell[4][i], iadd+my_union[iunion_scan].t.database.nvtxcell[5][i], iadd+my_union[iunion_scan].t.database.nvtxcell[6][i], iadd+my_union[iunion_scan].t.database.nvtxcell[7][i]);
#else
								// ������������ ������� ���� � ��������.
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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
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

	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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

// 10 ������ 2016. �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
// � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2. ��� ��������. ���� inumbercadr==0 �� ������ ����.
void exporttecplotxy360T_3D_part2_ianimation_series( integer maxelm, integer ncell, FLOW* &f, 
	TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, integer ikey,
	integer inumbercadr, doublereal time, BLOCK* b)
{
	const bool lite_export = true;
	// 16 ������ ����� ������� ��������� ������ �������,
	// ������ ���������� ����� ������.
	
	switch (ionly_solid_visible) {
	case WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE :
		std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE" << std::endl;
		break;
	case WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE :
		std::cout << "ionly_solid_visible =" << "WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE" << std::endl;
		break;
	}
	

	if (lite_export) {
		printf("lite export.\n");
	}
	else {
		printf("full export.\n");
	}

	// ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // ����� 1 ��� 3
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	doublereal* temp_shadow = nullptr;
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
		temp_shadow = new doublereal[t.maxelm + t.maxbound];
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			temp_shadow[i_1] = t.potent[i_1];
		}
	}

	doublereal** total_deformation_shadow = nullptr;
	if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
		total_deformation_shadow = new doublereal*[4];
		for (integer j_1 = 0; j_1 < 4; j_1++) {
			total_deformation_shadow[j_1] = new doublereal[t.maxelm + t.maxbound];
			for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
				total_deformation_shadow[j_1][i_1] = t.total_deformation[j_1][i_1];
			}
		}
	}


	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
#ifdef MINGW_COMPILLER
	int err = 0;
	if (inumbercadr == 0) {
		switch (ikey) {
		case 0: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		case 1: fp = fopen64("ALICEFLOW0_27_temp_animation.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		}
    }
	else {
		switch (ikey) {
		case 0: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		case 1: fp = fopen64("ALICEFLOW0_27_temp_animation.PLT", "ab");  break; // �� ��� ����� ��� ������ �������.
		default: fp = fopen64("ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		}
	}
	if (fp==NULL) err = 1;
#else
	errno_t err;
	if (inumbercadr == 0) {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_animation.PLT", "wb");  break; // �� ��� ����� ��� ������ �������.
		default: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "wb");  break;
		}
	}
	else {
		switch (ikey) {
		case 0: err = fopen_s(&fp, "ALICEFLOW0_07_temp_animation.PLT", "ab");  break;
		case 1: err = fopen_s(&fp, "ALICEFLOW0_27_temp_animation.PLT", "ab");  break; // �� ��� ����� ��� ������ �������.
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

		int c; // �������� ������
		integer ivarexport = 1; // �� ��������� ������ ���� ����������:
		

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


			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count > 0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (integer i = 0; i < flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // ������� ��� ����������� ������ �������
				}
			}

			if (ivarexport == 1) {

				integer ncell_shadow = 0;
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
					// ������������� ����������������� ��� ������� whot_is_block �����������
					// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
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
					// ������ ���������
					fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

					// ������ ��� ����������
					//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
					fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z, mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");  // ���������� ������ ���� ����������
					

				}
#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm, ncell_shadow);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n",time,   maxelm, ncell_shadow);
				}
#endif

				// ���������� ������� �����.
				//fprintf(fp, "$!ATTACHTEXT\n");
				//fprintf(fp, "POSITIONCOORDSYS = GRID\n ANCHORPOS\n {\n X=%e\n Y=%e\ Z=%e\n }\n",b[0].g.xS+0.05*(b[0].g.xE-b[0].g.xS), b[0].g.yS + 0.05*(b[0].g.yE - b[0].g.yS), b[0].g.zS + 0.05*(b[0].g.zE - b[0].g.zS));
				//fprintf(fp, "TEXTSHAPE\n {\n ISBOLD = NO\n SIZEUNITS = GRID\n HEIGHT = %e\n }\n",0.07*sqrt((b[0].g.xE - b[0].g.xS)*(b[0].g.xE - b[0].g.xS)+(b[0].g.yE - b[0].g.yS)*(b[0].g.yE - b[0].g.yS)+(b[0].g.zE - b[0].g.zS)*(b[0].g.zE - b[0].g.zS)));
				//fprintf(fp, "BOX\n{\nBOXTYPE = FILLED\n	MARGIN = 30\n COLOR = WHITE\n }\n");
				//fprintf(fp, "TEXT='time=%e s'\n\n", time);



				if (bvery_big_memory) {
					// extended print integer �� �������������.

					// ������ x
					for (integer i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (integer i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
					for (integer i = 0; i < t.database.maxelm; i++) {
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
					// ������ ���������
					fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");
				}

				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ������ �� ������
					ncell_shadow = 0;
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
					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
					{
						ncell_shadow = 0;
						for (integer i = 0; i < t.database.ncell; i++) {

							integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
							inode1 = t.database.nvtxcell[0][i] - 1;
							inode2 = t.database.nvtxcell[1][i] - 1;
							inode3 = t.database.nvtxcell[2][i] - 1;
							inode4 = t.database.nvtxcell[3][i] - 1;
							inode5 = t.database.nvtxcell[4][i] - 1;
							inode6 = t.database.nvtxcell[5][i] - 1;
							inode7 = t.database.nvtxcell[6][i] - 1;
							inode8 = t.database.nvtxcell[7][i] - 1;

							integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
							integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
							integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
							integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


							integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
							integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
							integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
							integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];


							integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
							integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
							integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
							integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];

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

							// ������������ ������ ������� ����.
							// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
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
							// � ��� ������ ���� ������� ����, ������� �� ���������� ������  ��������.
							ncell_shadow = ncell;
							ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;
						}
					}
				}
				if (inumbercadr == 0) {
					// ������ ����� ������� ������� � ������������� � �������������:
					
					//fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, log10_heat_flux_x, log10_heat_flux_y, log10_heat_flux_z,  mag_heat_flux, log10_mag_heat_flux, total_deformation, x_deformation, y_deformation, z_deformation\n");
						
					
				}
#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm, ncell_shadow);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"time=%e s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", time,  maxelm, ncell_shadow);
				}
#endif



				
					if (bvery_big_memory) {
						if (lite_export) {
							// ������ x
							for (integer i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// ������ y
							for (integer i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// ������ z
							for (integer i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.6f ", t.database.z[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
						}
						else {
							// ������ x
							for (integer i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", t.database.x[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// ������ y
							for (integer i = 0; i < t.database.maxelm; i++) {
								fprintf(fp, "%+.16f ", t.database.y[i]);
								if (i % 10 == 0) fprintf(fp, "\n");
							}
							// ������ z
							for (integer i = 0; i < t.database.maxelm; i++) {
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
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		

		// ������ ������ �����

		// ������ ���� ���������� ������������ ������.


		// ������ �����������
		
			if (lite_export) {
				for (integer i = 0; i < maxelm; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
						fprintf(fp, "%+.3f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.3f ", t.potent[i]);
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (integer i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
				for (integer i = 0; i < maxelm; i++) {
					if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
						fprintf(fp, "%+.16f ", temp_shadow[i]);
					}
					else {
						fprintf(fp, "%+.16f ", t.potent[i]);
					}
					if ( i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (integer i = maxelm; i < maxelm + t.maxbound; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
				for (integer i = 0; i < maxelm; i++) {
					fprintf(fp, "%+.4f ", t.prop[LAM][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (integer i = 0; i < t.maxbound; i++) {
						fprintf(fp, "%+.4f ", t.prop_b[LAM][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}
			else {
				for (integer i = 0; i < maxelm; i++) {
					fprintf(fp, "%+.16f ", t.prop[LAM][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					for (integer i = 0; i < t.maxbound; i++) {
						fprintf(fp, "%+.16f ", t.prop_b[LAM][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}
			}


			fprintf(fp, "\n");
		

		
			// ������ ����������������� ������� ���� ����������:
			if (ivarexport == 3) {
				// Speed
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						doublereal svx = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
						doublereal svy = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
						doublereal svz = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}


				if (bextendedprint) {
					// Speed
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						doublereal svx = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
						doublereal svy = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
						doublereal svz = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
						fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}


				fprintf(fp, "\n");

				// Pressure
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PRESS][t.ptr[0][i]]); // PRESSURE
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Pressure
					integer idfluid = 0;
					for (integer  i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[PRESS][i + maxelm]); // PRESSURE
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// PAM
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[PAM][t.ptr[0][i]]); // PAM
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// PAM
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[PAM][i + maxelm]); // PRESSURE
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// VX
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VX
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm]); // VX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// VY
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VY
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm]); // VY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// VZ
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// VZ
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm]); // VZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");


				// Rho
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
					}
					else fprintf(fp, "%+.16f ", t.prop[RHO][i]);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Rho
					for (integer i = 0; i < f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[0].prop_b[RHO][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// Mu
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Mu
					for (integer i = 0; i < f[0].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// Mut // ������������ ������������ ��������
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]] / f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);

						//fprintf(fp, "%+.16f ", 1.0*f[t.ptr[1][i]].icolor_different_fluid_domain[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// MUT
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						//fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm]); // MUT
						fprintf(fp, "%+.16f ", f[idfluid].potent[MUT][i + maxelm] / f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// ��� ������� ������ ������������� ��������� ����������� �������.
				// ��� ������������� ���������� �� ������ Distance_Wall.
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						//fprintf(fp, "%+.16f ", doublereal(i));
						if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
							(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY) ||
							(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
							(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
							(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
							(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						if ((f[0].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
							(f[0].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
							(f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
							(f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
							(f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
							(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
							fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
						}
						else fprintf(fp, "%+.16f ", 0.0);
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");



				// Curl // ������������ - ������ ������ ��������.
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[CURL][t.ptr[0][i]]); // CURL FBUF
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// Curl
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[CURL][i + maxelm]); // Curl
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// ������� ����������� �� ��������� �������� !!!.

				// GRADXVX
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVX
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVX][i + maxelm]); // GRADXVX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADYVX
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVX
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVX][i + maxelm]); // GRADYVX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADZVX
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVX][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVX
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVX][i + maxelm]); // GRADZVX
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADXVY
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVY
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVY][i + maxelm]); // GRADXVY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADYVY
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVY
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVY][i + maxelm]); // GRADYVY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADZVY
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVY][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVY
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADZVY][i + maxelm]); // GRADZVY
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADXVZ
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADXVZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADXVZ
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADXVZ][i + maxelm]); // GRADXVZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADYVZ
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADYVZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADYVZ
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
						fprintf(fp, "%+.16f ", f[idfluid].potent[GRADYVZ][i + maxelm]); // GRADYVZ
						if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
					}
				}

				fprintf(fp, "\n");

				// GRADZVZ
				for (integer i = 0; i < maxelm; i++) {
					if (t.ptr[1][i] > -1) {
						fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[GRADZVZ][t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				if (bextendedprint) {
					// GRADZVZ
					integer idfluid = 0;
					for (integer i = 0; i < f[idfluid].maxbound; i++) {
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

		// ������������� ����.
		for (integer i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// ���������� ����������.
		for (integer i = 0; i<t.maxelm; i++) {
			// ������ ���������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
		}

		for (integer i = 0; i<t.maxelm; i++) {
			// ������ ��������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
		}


		// ���������� � ����.
		
			if (lite_export) {
				doublereal buf0 = 0.0, buf1 = 0.0, buf2 = 0.0, buf3 = 0.0, buf4 = 0.0, buf5 = 0.0, buf6 = 0.0, buf7 = 0.0, buf8 = 0.0, buf9 = 0.0;

				// Heat Flux X
				for (integer i1 = 0; i1 < maxelm; i1++) {
					if ((i1 + 10) < maxelm) {

						integer i = i1;
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
						integer i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tx[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							integer i = i1;
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
							integer i = i1;
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

						integer i = i1;
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
						integer i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Ty[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							integer i = i1;
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
							integer i = i1;
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

						integer i = i1;
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
						integer i = i1;
						buf0 = signlog10(-t.prop[LAM][i] * Tz[i]);
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							integer i = i1;
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
							integer i = i1;
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
						integer i = i1;
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
						integer i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						fprintf(fp, "%+.6f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							integer i = i1;
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
						integer i = i1;
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
						integer i = i1;
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
							integer i = i1;
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

						integer i = i1;
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
						integer i = i1;
						buf0 = -t.prop[LAM][i] * Tx[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							integer i = i1;
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
							integer i = i1;
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

						integer i = i1;
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
						integer i = i1;
						buf0 = -t.prop[LAM][i] * Ty[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							integer i = i1;
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
							integer i = i1;
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

						integer i = i1;
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
						integer i = i1;
						buf0 = -t.prop[LAM][i] * Tz[i];
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {

							integer i = i1;
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
							integer i = i1;
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
						integer i = i1;
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
						integer i = i1;
						buf0 = sqrt((-t.prop[LAM][i] * Tx[i])*(-t.prop[LAM][i] * Tx[i]) + (-t.prop[LAM][i] * Ty[i])*(-t.prop[LAM][i] * Ty[i]) + (-t.prop[LAM][i] * Tz[i])*(-t.prop[LAM][i] * Tz[i]));
						fprintf(fp, "%+.16f ", buf0);
						if (i1 % 10 == 0) fprintf(fp, "\n");
					}
				}

				if (bextendedprint) {
					for (integer i1 = 0; i1 < t.maxbound; i1++) {
						if ((i1 + 10) < t.maxbound) {
							integer i = i1;
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
						integer i = i1;
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
						integer i = i1;
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
							integer i = i1;
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


		// ������������ ����������� ������.
		delete[] Tx;
		Tx = nullptr;

		delete[] Ty;
		Ty = nullptr;

		delete[] Tz;
		Tz = nullptr;

		
			fprintf(fp, "\n");

			for (integer j_6 = 0; j_6 < 4; j_6++) {

				// ������ ������ ����������
				if (lite_export) {
					for (integer i = 0; i < maxelm; i++) {


						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
						for (integer i = maxelm; i < maxelm + t.maxbound; i++) {
							if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
					for (integer i = 0; i < maxelm; i++) {
						if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
							fprintf(fp, "%e ", total_deformation_shadow[j_6][i]);
						}
						else {
							fprintf(fp, "%e ", t.total_deformation[j_6][i]);
						}
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (bextendedprint) {
						for (integer i = maxelm; i < maxelm + t.maxbound; i++) {
							if (ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) {
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
			// ������ ���������� � ���������� �����
			for (integer i = 0; i < t.database.ncell; i++) {
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
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

						integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
						integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
						integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
						integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


						integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
						integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
						integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
						integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



						integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
						integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
						integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
						integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];


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

						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
							// ������������ ������� ���� � ��������.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// ������������ ������� ���� � ��������.
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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
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

// ���������� ��� �������� ��� �������������� ������������� �������.
// 10 ������ 2016 . �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
// � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2.
void exporttecplotxy360T_3D_part2amg(TEMPER &t, doublereal* u, bool bextendedprint, integer imove)
{
	integer ianimate = 0;
	const integer flow_interior_count = 1;
	// imove 0 ��� 1 ��� ��������� � ���� ��� � �������.

	// ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // ����� 1 ��� 3
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp=fopen64("ALICEFLOW0_07_temp.PLT", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
#endif
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.


	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
	if ((err) != 0) {
		printf("Create File temp Error in function exporttecplotxy360T_3D_part2amg in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // �������� ������
		integer ivarexport = 1; // �� ��������� ������ ���� ����������:
		integer i = 0; // ������� �����

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


			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			
			// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
			for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
				ivarexport = 3; // ������� ��� ����������� ������ �������
			}
			

			if (ivarexport == 1) {
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// ������ ��� ����������
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
				fprintf(fp, "VARIABLES = x, y, z, u\n");  // ���������� ������ ������������ ����

				// ������ ���������� � �����
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
					// extended print integer �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
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
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				/*
				integer ncell_shadow = ncell;
				if (bsolid_static_only) {
					// ������ �� ������
					ncell_shadow = ncell;
				}
				else {
					if ((ionly_solid_visible == ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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
							// ������������ ������ ������� ����.
							// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
							if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
								&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {
								ncell_shadow++;
							}
						}
						if (ncell_shadow == 0) {
							// � ��� ������ ���� ������� ����, ������� �� ���������� ������  ��������.
							ncell_shadow = ncell;
							ionly_solid_visible = FLUID_AND_SOLID_BODY_VISIBLE;
						}
					}
				}
				*/
				// ������ ����� ������� ������� � ������������� � �������������:
				
				fprintf(fp, "\nVARIABLES = x, y, z, u\n");
				

				// ������ ���������� � �����
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
					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
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
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// ������ ������ �����

		// ������ ���� ���������� ������������ ������.

		/*
		// ������ �����������
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport == 3) {
			/*
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
					doublereal svy = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
					doublereal svz = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
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

			// Mut // ������������ ������������ ��������
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

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
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



			// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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

		// ������������� ����.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// ���������� ����������.
		for (i = 0; i<t.maxelm; i++) {
			// ������ ���������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, Tx, Ty, Tz);
		}

		for (i = 0; i<t.maxelm; i++) {
			// ������ ��������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, Tx, Ty, Tz);
		}


		// ���������� � ����.

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


		// ������������ ����������� ������.
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
			// ������ ���������� � ���������� �����
			for (i = 0; i < t.database.ncell; i++) {
				if (bsolid_static_only) {

#if doubleintprecision == 1
					// ������������ ������� ���� � ������ ��� � ������.
					fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
					// ������������ ������� ���� � ������ ��� � ������.
					fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						}
				else {
					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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
						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {

#if doubleintprecision == 1
							fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

								}
					}
					else {

#if doubleintprecision == 1
						// ������������ ������� ���� � ��������.
						fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// ������������ ������� ���� � ��������.
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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
	}
	//getchar();
	system("PAUSE");
	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

} // for amg


// 10 ������ 2016 . �������: ���� ������� ������ ������� ��������� �����, ����� �� ������� ���������� tecplot��,
// � �� ��� ������ � ��������� ������ ����� �������� ����� tecplot�� ���������� �� �������� ����������. 
// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2.
void exporttecplotxy360T_3D_part2_rev(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint, BLOCK* &b, integer lb)
{
	// ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage = false;

	FILE *fp=NULL;
	FILE *fp1=NULL; // ����� 1 ��� 3
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("ALICEFLOW0_07_temp.PLT", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
#endif
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.


	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	// w -write, b - binary.
	if ((err) != 0) {
		printf("Create File temp Error in exporttecplotxy360T_3D_part2_rev in my_export_tecplot3.c\n");
		//getchar();
		system("pause");

	}
	else {

		int c; // �������� ������
		integer ivarexport = 1; // �� ��������� ������ ���� ����������:
		integer i = 0; // ������� �����

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


			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count>0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (i = 0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport = 3; // ������� ��� ����������� ������ �������
				}
			}

			if (ivarexport == 1) {
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				// ������ ��� ����������
				//fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
				fprintf(fp, "VARIABLES = x, y, z, Temp, Lam, heat_flux_x, heat_flux_y, heat_flux_z, mag_heat_flux\n");  // ���������� ������ ���� ����������


#if doubleintprecision == 1		// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
			    }
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#else							// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}
#endif

				

				if (bvery_big_memory) {
					// extended print integer �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
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
				// ������ ���������
				fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

				integer ncell_shadow = ncell;
				if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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
						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1)
						&& (t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {
						    ncell_shadow++;
						}
						/*
						// ����� � �������������� �������������� ���.
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
						// � ��� ������ ���� ������� ����, ������� �� ���������� ������  ��������.
						ncell_shadow = ncell;
						ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;
					}
				}
				// ������ ����� ������� ������� � ������������� � �������������:
				fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z,  mag_heat_flux\n");
				

#if doubleintprecision == 1
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
			    }
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#else
				// ������ ���������� � �����
				if (bextendedprint) {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm + t.maxbound, ncell_shadow);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell_shadow);
				}
#endif

				
				if (bvery_big_memory) {
					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
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
				fclose(fp1); // �������� �����
			}
			if (bprintmessage) {
				printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}

		// ������ ������ �����

		// ������ ���� ���������� ������������ ������.


		// ������ �����������
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport == 3) {
			// Speed
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					doublereal svx = f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					doublereal svy = f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					doublereal svz = f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]] * f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
					fprintf(fp, "%+.16f ", sqrt(svx + svy + svz));
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}


			if (bextendedprint) {
				// Speed
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					doublereal svx = f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm];
					doublereal svy = f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm];
					doublereal svz = f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm] * f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm];
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VX
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_X_COMPONENT][i + maxelm]); // VX
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VY
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VY
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Y_COMPONENT][i + maxelm]); // VY
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// VZ
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// VZ
				integer idfluid = 0;
				for (i = 0; i < f[idfluid].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[idfluid].potent[VELOCITY_Z_COMPONENT][i + maxelm]); // VZ
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");


			// Rho
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[RHO][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					fprintf(fp, "%+.16f ", f[t.ptr[1][i]].prop[MU_DYNAMIC_VISCOSITY][t.ptr[0][i]]);
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
				}
				else fprintf(fp, "%+.16f ", 0.0);
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				// Mu
				for (i = 0; i < f[0].maxbound; i++) {
					fprintf(fp, "%+.16f ", f[0].prop_b[MU_DYNAMIC_VISCOSITY][i]);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			// Mut // ������������ ������������ ��������
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

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
			for (i = 0; i < maxelm; i++) {
				if (t.ptr[1][i] > -1) {
					//fprintf(fp, "%+.16f ", doublereal(i));
					if ((f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::SMAGORINSKY) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[t.ptr[1][i]].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
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
					if ((f[0].iflowregime == VISCOSITY_MODEL::ZEROEQMOD) ||
						(f[0].iflowregime == VISCOSITY_MODEL::SMAGORINSKY)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST) ||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS)||
						(f[0].iflowregime == VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST)) {
						fprintf(fp, "%+.16f ", f[idfluid].rdistWall[i + maxelm]); // Distance_Wall
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if ((i + maxelm) % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");



			// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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

		// ������������� ����.
		for (i = 0; i<t.maxelm + t.maxbound; i++) {
			Tx[i] = 0.0;
			Ty[i] = 0.0;
			Tz[i] = 0.0;
		}

		// ���������� ����������.
		for (i = 0; i<t.maxelm; i++) {
			// ������ ���������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
		}

		for (i = 0; i<t.maxelm; i++) {
			// ������ ��������� ����.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
		}


		// ���������� � ����.

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


		// ������������ ����������� ������.
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
			// ������ ���������� � ���������� �����
			for (i = 0; i < t.database.ncell; i++) {
				if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior>0))
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
					// ������������ ������ ������� ����.
					// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
					if ((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
						(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)) {

#if doubleintprecision == 1
						fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif

						}
				}
				else {


#if doubleintprecision == 1
					// ������������ ������� ���� � ��������.
					fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
					fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
					// ������������ ������� ���� � ��������.
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
					// ����������� ������� ����� � �������� ����
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // �������� �����
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}
		}

		// ���� �� ��� ������.
		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	//WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}//30 april 2016





// ��� ������ VARIATION Plot �� ANSYS Icepak.
// � 16 ������ 2018 ���� �� ������� �� ���������� ����Mesh_v0.42 ���������� �
// ����� �� ������� ����� ���������� ����������� (��������� ������� ����� �����
// ������� ����� ����� �������� � ����������� �����, ������� ������ ��������� � 
// ������������ ����� �� ���� ���������� ������������� ������� ���������.
void xyplot( FLOW* &fglobal, integer flow_interior, TEMPER &t) {
	
	
	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("xyplot.txt", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "xyplot.txt", "w");
#endif
	

	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else if (fglobal[0].maxelm>0) {

		if (fp != NULL) {

			
			doublereal x = 0.0e-3, y = 0.0e-3, z = 0.0e-3; // ����� ����� ������� �������� �����
			x = Tochka_position_X0_for_XY_Plot;
			y = Tochka_position_Y0_for_XY_Plot;
			z = Tochka_position_Z0_for_XY_Plot;

			const integer ifi = 0;

			TOCHKA sp0, sp1; // ������� [sp0, sp1].
			doublereal xmin = 1.0e30, xmax = -1.0e30;
			doublereal ymin = 1.0e30, ymax = -1.0e30;
			doublereal zmin = 1.0e30, zmax = -1.0e30;
			for (integer iP = 0; iP < fglobal[ifi].maxelm; iP++) {
				TOCHKA p;
				p.x = 0;
				p.y = 0;
				p.z = 0;
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100); // ���������� ��������� ������ ��.
				if (p.x < xmin) xmin = p.x;
				if (p.x > xmax) xmax = p.x;
				if (p.y < ymin) ymin = p.y;
				if (p.y > ymax) ymax = p.y;
				if (p.z < zmin) zmin = p.z;
				if (p.z > zmax) zmax = p.z;
			}

			
			integer iplane = XZ_PLANE; // ��������� ���������������� �����.
			switch (idirectional_for_XY_Plot) {
			case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // X
				iplane = YZ_PLANE; // ��������� ���������������� �����.
				sp0.y = sp1.y = y;
				sp0.z = sp1.z = z;
				sp0.x = xmin;
				sp1.x = xmax;
				break;
			case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: // Y
				iplane = XZ_PLANE; // ��������� ���������������� �����.
				sp0.x = sp1.x = x;
				sp0.z = sp1.z = z;
				sp0.y = ymin;
				sp1.y = ymax;
				break;
			case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // Z
				iplane = XY_PLANE; // ��������� ���������������� �����.
				sp0.y = sp1.y = y;
				sp0.x = sp1.x = x;
				sp0.z = zmin;
				sp1.z = zmax;
				break;
			default: // X
				iplane = YZ_PLANE; // ��������� ���������������� �����.
				sp0.y = sp1.y = y;
				sp0.z = sp1.z = z;
				sp0.x = xmin;
				sp1.x = xmax;
				break;
			}


			integer G=E_SIDE;
			switch (iplane) {
			case XY_PLANE: G = T_SIDE;  break;
			case XZ_PLANE: G = N_SIDE;  break;
			case YZ_PLANE: G = E_SIDE;  break;
			}

			doublereal dopusk = 1.0e30;
			for (integer iP = 0; iP <  fglobal[ifi].maxelm; iP++) {
				TOCHKA p;
				p.x = 0;
				p.y = 0;
				p.z = 0;
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100); // ���������� ��������� ������ ��.
				if (distance_Point_to_Segment2(p, sp0, sp1) < kR_size_cell * volume3D_ray_tracing(iP, fglobal[ifi].nvtx, fglobal[ifi].pa)) {
					// ����� ����������� �������.
					if (kR_size_cell * volume3D_ray_tracing(iP, t.nvtx, t.pa) < dopusk) {
						dopusk = 0.49 * volume3D_ray_tracing(iP, t.nvtx, t.pa);
					}

				}
			}

			SORT_XY_PLOT_FLOW* xy_arr = new SORT_XY_PLOT_FLOW[fglobal[ifi].maxelm];

			integer icounter = -1;

			for (integer iP = 0; iP < fglobal[ifi].maxelm; iP++) {
				TOCHKA p;
				p.x = 0;
				p.y = 0;
				p.z = 0;
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100); // ���������� ��������� ������ ��.

				if (distance_Point_to_Segment2(p, sp0, sp1) < dopusk) {
					// ����� ����������� �������.

					icounter++;
					xy_arr[icounter].Vx = fglobal[ifi].potent[VELOCITY_X_COMPONENT][iP];
					xy_arr[icounter].Vy = fglobal[ifi].potent[VELOCITY_Y_COMPONENT][iP];
					xy_arr[icounter].Vz = fglobal[ifi].potent[VELOCITY_Z_COMPONENT][iP];
					xy_arr[icounter].Speed = sqrt((fglobal[ifi].potent[VELOCITY_X_COMPONENT][iP])*(fglobal[ifi].potent[VELOCITY_X_COMPONENT][iP])+
						(fglobal[ifi].potent[VELOCITY_Y_COMPONENT][iP])*(fglobal[ifi].potent[VELOCITY_Y_COMPONENT][iP])+
						(fglobal[ifi].potent[VELOCITY_Z_COMPONENT][iP])*(fglobal[ifi].potent[VELOCITY_Z_COMPONENT][iP]));
					xy_arr[icounter].Pam = fglobal[ifi].potent[PAM][iP];
					xy_arr[icounter].Press = fglobal[ifi].potent[PRESS][iP];
					xy_arr[icounter].Fbuf = fglobal[ifi].potent[FBUF][iP];
					xy_arr[icounter].massflux_in_gran = fglobal[ifi].mf[iP][G];
					// ����������� �� ��������.
					xy_arr[icounter].GradXVx = fglobal[ifi].potent[GRADXVX][iP];
					xy_arr[icounter].GradYVx = fglobal[ifi].potent[GRADYVX][iP];
					xy_arr[icounter].GradZVx = fglobal[ifi].potent[GRADZVX][iP];
					xy_arr[icounter].GradXVy = fglobal[ifi].potent[GRADXVY][iP];
					xy_arr[icounter].GradYVy = fglobal[ifi].potent[GRADYVY][iP];
					xy_arr[icounter].GradZVy = fglobal[ifi].potent[GRADZVY][iP];
					xy_arr[icounter].GradXVz = fglobal[ifi].potent[GRADXVZ][iP];
					xy_arr[icounter].GradYVz = fglobal[ifi].potent[GRADYVZ][iP];
					xy_arr[icounter].GradZVz = fglobal[ifi].potent[GRADZVZ][iP];
					//***
					xy_arr[icounter].Curl = fglobal[ifi].potent[CURL][iP];
					xy_arr[icounter].GradXPress = fglobal[ifi].potent[GRADXPRESS][iP];
					xy_arr[icounter].GradYPress = fglobal[ifi].potent[GRADYPRESS][iP];
					xy_arr[icounter].GradZPress = fglobal[ifi].potent[GRADZPRESS][iP];
					xy_arr[icounter].Nusha = fglobal[ifi].potent[NUSHA][iP];
					xy_arr[icounter].GradXNusha = fglobal[ifi].potent[GRADXNUSHA][iP];
					xy_arr[icounter].GradYNusha = fglobal[ifi].potent[GRADYNUSHA][iP];
					xy_arr[icounter].GradZNusha = fglobal[ifi].potent[GRADZNUSHA][iP];
					//++++
					xy_arr[icounter].Ke = fglobal[ifi].potent[TURBULENT_KINETIK_ENERGY][iP];
					xy_arr[icounter].Omega = fglobal[ifi].potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
					xy_arr[icounter].GradXKe = fglobal[ifi].potent[GRADXTURBULENT_KINETIK_ENERGY][iP];
					xy_arr[icounter].GradYKe = fglobal[ifi].potent[GRADYTURBULENT_KINETIK_ENERGY][iP];
					xy_arr[icounter].GradZKe = fglobal[ifi].potent[GRADZTURBULENT_KINETIK_ENERGY][iP];
					xy_arr[icounter].GradXOmega = fglobal[ifi].potent[GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
					xy_arr[icounter].GradYOmega = fglobal[ifi].potent[GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];
					xy_arr[icounter].GradZOmega = fglobal[ifi].potent[GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][iP];


					switch (idirectional_for_XY_Plot) {
					case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
						xy_arr[icounter].x_argument = p.x;
						xy_arr[icounter].GradPam = fglobal[ifi].potent[GRADZPAM][iP];
						xy_arr[icounter].GradPress = fglobal[ifi].potent[GRADZPRESS][iP];						
						break;
					case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: //XZ
						xy_arr[icounter].x_argument = p.y;
						xy_arr[icounter].GradPam = fglobal[ifi].potent[GRADYPAM][iP];
						xy_arr[icounter].GradPress = fglobal[ifi].potent[GRADYPRESS][iP];
						break;
					case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
						xy_arr[icounter].x_argument = p.z;
						xy_arr[icounter].GradPam = fglobal[ifi].potent[GRADZPAM][iP];
						xy_arr[icounter].GradPress = fglobal[ifi].potent[GRADZPRESS][iP];
						break;
					default:
						// YZ
						xy_arr[icounter].x_argument = p.x;
						xy_arr[icounter].GradPam = fglobal[ifi].potent[GRADZPAM][iP];
						xy_arr[icounter].GradPress = fglobal[ifi].potent[GRADZPRESS][iP];
						break;
					}
				}
			}

			// ���������� �� �����������.
			std::sort(xy_arr, xy_arr + icounter + 1, compare_XY_PLOT_FLOW);		
			
			// ���� ��������� �������� � ����� � ��� �� 
			// ����������, �� �� ��������� ��������� �������,
			// �������� ������ � ����� ������ ��������.
			doublereal pos = xy_arr[0].x_argument;
			integer jposY = 0, jmaximumPos = jposY;
			
			doublereal maximumTemp = xy_arr[jposY].Speed;// ������ ��������.
			if (0) {
				for (integer j = 1; j <= icounter; j++) {
					if (pos == xy_arr[j].x_argument) {
						xy_arr[j] = xy_arr[jposY];
					}
					else {
						jposY = j;
						pos = xy_arr[j].x_argument;
					}
				}
			}
			else {
				for (integer j = 1; j <= icounter; j++) {
					if (pos == xy_arr[j].x_argument) {
						if (xy_arr[j].Speed > maximumTemp) {
							maximumTemp = xy_arr[j].Speed;
							jmaximumPos = j;
						}
					}
					else {
						// ������ ����������� �������� �����������
						// � �����������, ��� ���� �������������.
						for (integer j31 = jposY; j31 < j; j31++) {
							xy_arr[j31] = xy_arr[jmaximumPos];
						}
						jposY = j;
						maximumTemp = xy_arr[jposY].Speed;// ������ ��������.
						jmaximumPos = jposY;
						pos = xy_arr[j].x_argument;
					}
				}
			}

			fprintf(fp, "Xo=%e, Yo=%e, Zo=%e, ", x, y, z);
			switch (idirectional_for_XY_Plot) {
			case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
				fprintf(fp, "directional=X\n");
				break;
			case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: //XZ
				fprintf(fp, "directional=Y\n");
				break;
			case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
				fprintf(fp, "directional=Z\n");
				break;
			default:
				fprintf(fp, "directional=X\n");
				break;
			}

			fprintf(fp, "position,\tVx,\tVy,\tVz,\tSpeed,\tPam,\tPress,\tFbuf,\tGRADPRESS,\tGRADPAM,\tmassfluxingran");
			fprintf(fp, ", \tGradXVx, \tGradYVx, \tGradZVx, \tGradXVy, \tGradYVy, \tGradZVy, \tGradXVz, \tGradYVz, \tGradZVz");
			fprintf(fp, ", \tCurl, \tGradXPress,  \tGradYPress,  \tGradZPress, \tNusha_turbulence, \tGradXNusha, \tGradYNusha, \tGradZNusha");
			fprintf(fp, ", \tkinetik_energy, \tomega, \tGradXKe, \tGradYKe, \tGradZKe, \tGradXOmega, \tGradYOmega, \tGradZOmega\n");
			for (integer j = 0; j <= icounter; j++) {
				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f",
					xy_arr[j].x_argument, xy_arr[j].Vx, xy_arr[j].Vy, xy_arr[j].Vz,
					xy_arr[j].Speed, xy_arr[j].Pam, xy_arr[j].Press, xy_arr[j].Fbuf,
					xy_arr[j].GradPress, xy_arr[j].GradPam, xy_arr[j].massflux_in_gran);
				fprintf(fp," %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f",
					xy_arr[j].GradXVx, xy_arr[j].GradYVx, xy_arr[j].GradZVx,
					xy_arr[j].GradXVy, xy_arr[j].GradYVy, xy_arr[j].GradZVy,
					xy_arr[j].GradXVz, xy_arr[j].GradYVz, xy_arr[j].GradZVz);
				fprintf(fp, " %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f",
					xy_arr[j].Curl, xy_arr[j].GradXPress, xy_arr[j].GradYPress,
					xy_arr[j].GradZPress, xy_arr[j].Nusha, xy_arr[j].GradXNusha,
					xy_arr[j].GradYNusha, xy_arr[j].GradZNusha);
				fprintf(fp, " %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
					xy_arr[j].Ke, xy_arr[j].Omega, xy_arr[j].GradXKe,
					xy_arr[j].GradYKe, xy_arr[j].GradZKe, xy_arr[j].GradXOmega,
					xy_arr[j].GradYOmega, xy_arr[j].GradZOmega);

			}

			delete[] xy_arr;			
			
			fclose(fp);
		}
	}
} // xyplot()


// ��� ������ VARIATION Plot �� ANSYS Icepak.
// � 16 ������ 2018 ���� �� ������� �� ���������� ����Mesh_v0.42 ���������� �
// ����� �� ������� ����� ���������� ����������� (��������� ������� ����� �����
// ������� ����� ����� �������� � ����������� �����, ������� ������ ��������� � 
// ������������ ����� �� ���� ���������� ������������� ������� ���������.
void xyplot2(FLOW* &fglobal, integer flow_interior, TEMPER &t) {
	FILE *fp = NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("xyplot.txt", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "xyplot.txt", "w");
#endif


	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else {

		TOCHKA p;
		doublereal epsilon = 1e30;
		doublereal dist;
		doublereal x = 0.0e-3, y = 0.0e-3, z = 0.0e-3; // ����� ����� ������� �������� �����
		x = Tochka_position_X0_for_XY_Plot;
		y = Tochka_position_Y0_for_XY_Plot;
		z = Tochka_position_Z0_for_XY_Plot;

		integer ifi = 0, iPf = 0;
		integer iplane = XZ_PLANE; // ��������� ���������������� �����.
		switch (idirectional_for_XY_Plot) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // X
			iplane = YZ_PLANE; // ��������� ���������������� �����.
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: // Y
			iplane = XZ_PLANE; // ��������� ���������������� �����.
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // Z
			iplane = XY_PLANE; // ��������� ���������������� �����.
			break;
		default: // X
			iplane = YZ_PLANE; // ��������� ���������������� �����.
			break;
		}


		for (integer i = 0; i<flow_interior; i++) {
			for (integer iP = 0; iP<fglobal[i].maxelm; iP++) {
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100); // ���������� ��������� ������ ��.
				dist = sqrt(fabs(x - p.x)*fabs(x - p.x) + fabs(y - p.y)*fabs(y - p.y) + fabs(z - p.z)*fabs(z - p.z));
				if (dist<epsilon) {
					epsilon = dist;
					ifi = i;
					iPf = iP;
				}
			}
		}
		// ��������� � ������ 
		switch (iplane) {
		case XY_PLANE: while (fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]<fglobal[ifi].maxelm) iPf = fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]; break;
		case XZ_PLANE: while (fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]<fglobal[ifi].maxelm) iPf = fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]; break;
		case YZ_PLANE: while (fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]<fglobal[ifi].maxelm) iPf = fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]; break;
		}

		integer G;
		switch (iplane) {
		case XY_PLANE: G = T_SIDE;  break;
		case XZ_PLANE: G = N_SIDE;  break;
		case YZ_PLANE: G = E_SIDE;  break;
		}


		fprintf(fp, "position,\tVx,\tVy,\tVz,\tPam,\tPress,\tFbuf,\tGRADPRESS,\tGRADPAM,\tmassfluxingran\n");
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
		//volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
		volume3D_q(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
		center_cord3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100); // ���������� ��������� ������ ��.
		switch (iplane) {
		case XY_PLANE: fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			p.z - 0.5*dz, fglobal[ifi].potent[VELOCITY_X_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[VELOCITY_Y_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[VELOCITY_Z_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[GRADZPRESS][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].potent[GRADZPAM][fglobal[ifi].neighbors_for_the_internal_node[B_SIDE][0][iPf]],
			fglobal[ifi].mf[iPf][G]);
			break;
		case XZ_PLANE:  fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			p.y - 0.5*dy, fglobal[ifi].potent[VELOCITY_X_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[VELOCITY_Y_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[VELOCITY_Z_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[GRADYPRESS][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].potent[GRADYPAM][fglobal[ifi].neighbors_for_the_internal_node[S_SIDE][0][iPf]],
			fglobal[ifi].mf[iPf][G]);
			break;
		case YZ_PLANE:  fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
			p.x - 0.5*dx, fglobal[ifi].potent[VELOCITY_X_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[VELOCITY_Y_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[VELOCITY_Z_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[GRADXPRESS][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].potent[GRADXPAM][fglobal[ifi].neighbors_for_the_internal_node[W_SIDE][0][iPf]],
			fglobal[ifi].mf[iPf][G]);
			break;
		}
		switch (iplane) {
		case XY_PLANE: while (iPf<fglobal[ifi].maxelm) {
			center_cord3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100);
			fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
				p.z, fglobal[ifi].potent[VELOCITY_X_COMPONENT][iPf],
				fglobal[ifi].potent[VELOCITY_Y_COMPONENT][iPf],
				fglobal[ifi].potent[VELOCITY_Z_COMPONENT][iPf],
				fglobal[ifi].potent[PAM][iPf],
				fglobal[ifi].potent[PRESS][iPf],
				fglobal[ifi].potent[FBUF][iPf],
				fglobal[ifi].potent[GRADZPRESS][iPf],
				fglobal[ifi].potent[GRADZPAM][iPf],
				fglobal[ifi].mf[iPf][G]);
			if (fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf] >= fglobal[ifi].maxelm) {
				volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
					p.z + 0.5*dz, fglobal[ifi].potent[VELOCITY_X_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[VELOCITY_Y_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[VELOCITY_Z_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[GRADZPRESS][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].potent[GRADZPAM][fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf]],
					fglobal[ifi].mf[iPf][G]);
			}
			iPf = fglobal[ifi].neighbors_for_the_internal_node[T_SIDE][0][iPf];
		} break;
		case XZ_PLANE: while (iPf<fglobal[ifi].maxelm) {
			center_cord3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100);
			fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f  %+.16f %+.16f %+.16f\n",
				p.y, fglobal[ifi].potent[VELOCITY_X_COMPONENT][iPf],
				fglobal[ifi].potent[VELOCITY_Y_COMPONENT][iPf],
				fglobal[ifi].potent[VELOCITY_Z_COMPONENT][iPf],
				fglobal[ifi].potent[PAM][iPf],
				fglobal[ifi].potent[PRESS][iPf],
				fglobal[ifi].potent[FBUF][iPf],
				fglobal[ifi].potent[GRADYPRESS][iPf],
				fglobal[ifi].potent[GRADYPAM][iPf],
				fglobal[ifi].mf[iPf][G]);
			if (fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf] >= fglobal[ifi].maxelm) {
				volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
					p.y + 0.5*dy, fglobal[ifi].potent[VELOCITY_X_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[VELOCITY_Y_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[VELOCITY_Z_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[GRADYPRESS][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].potent[GRADYPAM][fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf]],
					fglobal[ifi].mf[iPf][G]);
			}
			/*

			#if doubleintprecision == 1
			// ����� ������������������ ����� ��� �������.
			printf("iPf=%lld\n",iPf);
			if (fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
			printf("iPffinish=%lld\n",fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1);
			getchar();
			}
			#else
			// ����� ������������������ ����� ��� �������.
			printf("iPf=%d\n",iPf);
			if (fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
			printf("iPffinish=%d\n",fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1);
			getchar();
			}
			#endif


			*/
			iPf = fglobal[ifi].neighbors_for_the_internal_node[N_SIDE][0][iPf];
		} break;
		case YZ_PLANE: while (iPf<fglobal[ifi].maxelm) {
			center_cord3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, p, 100);
			fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
				p.y, fglobal[ifi].potent[VELOCITY_X_COMPONENT][iPf],
				fglobal[ifi].potent[VELOCITY_Y_COMPONENT][iPf],
				fglobal[ifi].potent[VELOCITY_Z_COMPONENT][iPf],
				fglobal[ifi].potent[PAM][iPf],
				fglobal[ifi].potent[PRESS][iPf],
				fglobal[ifi].potent[FBUF][iPf],
				fglobal[ifi].potent[GRADXPRESS][iPf],
				fglobal[ifi].potent[GRADXPAM][iPf],
				fglobal[ifi].mf[iPf][G]);
			if (fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf] >= fglobal[ifi].maxelm) {
				volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
				fprintf(fp, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n",
					p.x + 0.5*dx, fglobal[ifi].potent[VELOCITY_X_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[VELOCITY_Y_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[VELOCITY_Z_COMPONENT][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[PAM][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[PRESS][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[FBUF][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[GRADXPRESS][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].potent[GRADXPAM][fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf]],
					fglobal[ifi].mf[iPf][G]);
			}
			iPf = fglobal[ifi].neighbors_for_the_internal_node[E_SIDE][0][iPf];
		} break;
		}
		fclose(fp);
	}
} // xyplot2 ���������� ������.

// ���������� ������� ����������� ����� ����� 
void xyplot_temp2(TEMPER &t, doublereal* tempfiltr) {

	// tempfiltr - ������������ ���������� ������������� �����������.
	// ��� ������������ ����� xyplotT.txt.

	FILE *fp=NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("xyplotT.txt", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "xyplotT.txt", "w");
#endif
	
	FILE* fp1 = NULL;

#ifdef MINGW_COMPILLER
	int err1 = 0;
	fp1 = fopen64("xyplotT1.PLT", "w");
	if (fp1 == NULL) err1 = 1;
#else
	errno_t err1 = 0;
	err1 = fopen_s(&fp1, "xyplotT1.PLT", "w");
#endif

	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else if (t.maxelm > 0) {
		if (fp != NULL)
		{

			// �������� ! ��������� ������� ����� ����� ������� ����� ��������� ����� � ��������� ������� ������ ����� ����� ���������������.
			TOCHKA p;
			p.x = 0;
			p.y = 0;
			p.z = 0;
			doublereal epsilon = 1e30;
			doublereal dist;
			//doublereal x = 2.0245e-3, y = 0.2268e-3, z = 0.0125e-3; // ����� ����� ������� �������� �����
			doublereal x = 3e-3, y = -0.00005e-3, z = 0.57425e-3;
			x = Tochka_position_X0_for_XY_Plot;
			y = Tochka_position_Y0_for_XY_Plot;
			z = Tochka_position_Z0_for_XY_Plot;
			integer iPf = 0;
			integer iplane = YZ_PLANE; // ��������� ���������������� �����.
			switch (idirectional_for_XY_Plot) {
			case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
				iplane = YZ_PLANE; // ��������� ���������������� �����.
				break;
			case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: //XZ
				iplane = XZ_PLANE; // ��������� ���������������� �����.
				break;
			case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
				iplane = XY_PLANE; // ��������� ���������������� �����.
				break;
			default:
				iplane = YZ_PLANE; // ��������� ���������������� �����.
				break;
			}

			for (integer iP = 0; iP < t.maxelm; iP++) {
				center_cord3D(iP, t.nvtx, t.pa, p,100); // ���������� ��������� ������ ��.
				dist = sqrt(fabs(x - p.x)*fabs(x - p.x) + fabs(y - p.y)*fabs(y - p.y) + fabs(z - p.z)*fabs(z - p.z));
				if (dist < epsilon) {
					epsilon = dist;
					iPf = iP;
				}
			}

			// ��������� � ������ 
			if (t.maxelm > 0) {
				switch (iplane) {
				case XY_PLANE: while (t.neighbors_for_the_internal_node[B_SIDE][0][iPf] < t.maxelm) iPf = t.neighbors_for_the_internal_node[B_SIDE][0][iPf]; break;
				case XZ_PLANE: while (t.neighbors_for_the_internal_node[S_SIDE][0][iPf] < t.maxelm) iPf = t.neighbors_for_the_internal_node[S_SIDE][0][iPf]; break;
				case YZ_PLANE: while (t.neighbors_for_the_internal_node[W_SIDE][0][iPf] < t.maxelm) iPf = t.neighbors_for_the_internal_node[W_SIDE][0][iPf]; break;
				}
			}

			//integer G;
			//switch (iplane) {
			//case XY_PLANE: G = T_SIDE;  break;
			//case XZ_PLANE: G = N_SIDE;  break;
			//case YZ_PLANE: G = E_SIDE;  break;
			//}


			fprintf(fp, "position,\ttemperature,\ttemperature_avg\n");
			fprintf(fp1, "VARIABLE= position temperature\n");
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
			if (t.maxelm > 0) {
				volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
				center_cord3D(iPf, t.nvtx, t.pa, p, 100); // ���������� ��������� ������ ��.
			}
			switch (iplane) {
			case XY_PLANE: fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.z - 0.5*dz, t.potent[t.neighbors_for_the_internal_node[B_SIDE][0][iPf]],
				tempfiltr[t.neighbors_for_the_internal_node[B_SIDE][0][iPf]]);

				fprintf(fp1, "%+.16f %+.16f\n",
					p.z - 0.5 * dz, t.potent[t.neighbors_for_the_internal_node[B_SIDE][0][iPf]]);
				break;
			case XZ_PLANE:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.y - 0.5*dy, t.potent[t.neighbors_for_the_internal_node[S_SIDE][0][iPf]],
				tempfiltr[t.neighbors_for_the_internal_node[S_SIDE][0][iPf]]);

				fprintf(fp1, "%+.16f %+.16f\n",
					p.y - 0.5 * dy, t.potent[t.neighbors_for_the_internal_node[S_SIDE][0][iPf]]);
				break;
			case YZ_PLANE:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.x - 0.5*dx, t.potent[t.neighbors_for_the_internal_node[W_SIDE][0][iPf]],
				tempfiltr[t.neighbors_for_the_internal_node[W_SIDE][0][iPf]]);

				fprintf(fp1, "%+.16f %+.16f\n",
					p.x - 0.5 * dx, t.potent[t.neighbors_for_the_internal_node[W_SIDE][0][iPf]]);
				break;
			}
			switch (iplane) {
			case XY_PLANE: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p,100);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.z, t.potent[iPf],
					tempfiltr[iPf]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.z, t.potent[iPf]);
				if (t.neighbors_for_the_internal_node[T_SIDE][0][iPf] >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.z + 0.5*dz, t.potent[t.neighbors_for_the_internal_node[T_SIDE][0][iPf]],
						tempfiltr[t.neighbors_for_the_internal_node[T_SIDE][0][iPf]]);
					fprintf(fp1, "%+.16f %+.16f\n",
						p.z + 0.5 * dz, t.potent[t.neighbors_for_the_internal_node[T_SIDE][0][iPf]]);

				}
				iPf = t.neighbors_for_the_internal_node[T_SIDE][0][iPf];
			} break;
			case XZ_PLANE: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p,100);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.y, t.potent[iPf],
					tempfiltr[iPf]);

				fprintf(fp1, "%+.16f %+.16f\n",
					p.y, t.potent[iPf]);

				if (t.neighbors_for_the_internal_node[N_SIDE][0][iPf] >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.y + 0.5*dy, t.potent[t.neighbors_for_the_internal_node[N_SIDE][0][iPf]],
						tempfiltr[t.neighbors_for_the_internal_node[N_SIDE][0][iPf]]);

					fprintf(fp1, "%+.16f %+.16f\n",
						p.y + 0.5 * dy, t.potent[t.neighbors_for_the_internal_node[N_SIDE][0][iPf]]);

				}
				/*

				#if doubleintprecision == 1
					// ����� ������������������ ����� ��� �������.
					printf("iPf=%lld\n",iPf);
					if (fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
						printf("iPffinish=%lld\n",fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1);
						getchar();
					}
				#else
					// ����� ������������������ ����� ��� �������.
					printf("iPf=%d\n",iPf);
					if (fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1>=fglobal[ifi].maxelm) {
						printf("iPffinish=%d\n",fglobal[ifi].neighbors_for_the_internal_node[NSIDE][iPf].iNODE1);
						getchar();
					}
				#endif

				
				*/
				iPf = t.neighbors_for_the_internal_node[N_SIDE][0][iPf];
			} break;
			case YZ_PLANE: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p,100);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.x, t.potent[iPf],
					tempfiltr[iPf]);

				fprintf(fp1, "%+.16f %+.16f\n",
					p.x, t.potent[iPf]);

				if (t.neighbors_for_the_internal_node[E_SIDE][0][iPf] >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.x + 0.5*dx, t.potent[t.neighbors_for_the_internal_node[E_SIDE][0][iPf]],
						tempfiltr[t.neighbors_for_the_internal_node[E_SIDE][0][iPf]]);

					fprintf(fp1, "%+.16f %+.16f\n",
						p.x + 0.5 * dx, t.potent[t.neighbors_for_the_internal_node[E_SIDE][0][iPf]]);
				}
				iPf = t.neighbors_for_the_internal_node[E_SIDE][0][iPf];
			} break;
			}
			fclose(fp);
			fclose(fp1);
		}
	}
} // xyplot_temp2



typedef struct TSORT_XY_PLOT {
	doublereal x_argument;
	doublereal y_function1;
	doublereal y_function2;
} SORT_XY_PLOT;

bool compare_XY_PLOT(SORT_XY_PLOT& Amat1, SORT_XY_PLOT& Amat2) {
	return (Amat1.x_argument < Amat2.x_argument);
} // compare_XY_PLOT

  // ���������� ������� ����������� ����� ����� 
// �������� � �� ���� �����. 18,05,2020.
void xyplot_temp(TEMPER &t, doublereal* tempfiltr) {

	// tempfiltr - ������������ ���������� ������������� �����������.
	// ��� ������������ ����� xyplotT.txt.

	FILE *fp = NULL;
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("xyplotT.txt", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "xyplotT.txt", "w");
#endif


	FILE* fp1 = NULL;

#ifdef MINGW_COMPILLER
	int err1 = 0;
	fp1 = fopen64("xyplotT1.PLT", "w");
	if (fp1 == NULL) err1 = 1;
#else
	errno_t err1 = 0;
	err1 = fopen_s(&fp1, "xyplotT1.PLT", "w");
#endif

	if ((err) != 0) {
		printf("Create File xyplot Error\n");
		//getchar();
		system("pause");

	}
	else if (t.maxelm > 0) {
		if (fp != NULL)
		{

			// �������� ! ��������� ������� ����� ����� ������� ����� ��������� ����� � ��������� ������� ������ ����� ����� ���������������.
			
			//doublereal x = 2.0245e-3, y = 0.2268e-3, z = 0.0125e-3; // ����� ����� ������� �������� �����
			doublereal x = 3e-3, y = -0.00005e-3, z = 0.57425e-3;
			x = Tochka_position_X0_for_XY_Plot;
			y = Tochka_position_Y0_for_XY_Plot;
			z = Tochka_position_Z0_for_XY_Plot;

			
			TOCHKA sp0, sp1; // ������� [sp0, sp1].
			doublereal xmin = 1.0e30, xmax = -1.0e30;
			doublereal ymin = 1.0e30, ymax = -1.0e30;
			doublereal zmin = 1.0e30, zmax = -1.0e30;
			for (integer iP = 0; iP < t.maxelm; iP++) {
				TOCHKA p;
				p.x = 0;
				p.y = 0;
				p.z = 0;
				center_cord3D(iP, t.nvtx, t.pa, p, 100); // ���������� ��������� ������ ��.
				if (p.x < xmin) xmin = p.x;
				if (p.x > xmax) xmax = p.x;
				if (p.y < ymin) ymin = p.y;
				if (p.y > ymax) ymax = p.y;
				if (p.z < zmin) zmin = p.z;
				if (p.z > zmax) zmax = p.z;
			}
			
			integer iplane = YZ_PLANE; // ��������� ���������������� �����.
			switch (idirectional_for_XY_Plot) {
			case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
				iplane = YZ_PLANE; // ��������� ���������������� �����.
				sp0.y = sp1.y = y;
				sp0.z = sp1.z = z;
				sp0.x = xmin;
				sp1.x = xmax;
				break;
			case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: //XZ
				iplane = XZ_PLANE; // ��������� ���������������� �����.
				sp0.x = sp1.x = x;
				sp0.z = sp1.z = z;
				sp0.y = ymin;
				sp1.y = ymax;
				break;
			case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
				iplane = XY_PLANE; // ��������� ���������������� �����.
				sp0.y = sp1.y = y;
				sp0.x = sp1.x = x;
				sp0.z = zmin;
				sp1.z = zmax;
				break;
			default:
				iplane = YZ_PLANE; // ��������� ���������������� �����.
				sp0.y = sp1.y = y;
				sp0.z = sp1.z = z;
				sp0.x = xmin;
				sp1.x = xmax;
				break;
			}

			
			doublereal dopusk = 1.0e30;
			for (integer iP = 0; iP < t.maxelm; iP++) {
				TOCHKA p;
				p.x = 0;
				p.y = 0;
				p.z = 0;
				center_cord3D(iP, t.nvtx, t.pa, p, 100); // ���������� ��������� ������ ��.
				if (distance_Point_to_Segment2(p, sp0, sp1) < kR_size_cell * volume3D_ray_tracing(iP, t.nvtx, t.pa)) {
					// ����� ����������� �������.
					if (kR_size_cell * volume3D_ray_tracing(iP, t.nvtx, t.pa) < dopusk) {
						dopusk = 0.49 * volume3D_ray_tracing(iP, t.nvtx, t.pa);
					}

				}
			}

			
 
			SORT_XY_PLOT* xy_arr = new SORT_XY_PLOT[t.maxelm];
			
			integer icounter = -1;

			for (integer iP = 0; iP < t.maxelm; iP++) {
				TOCHKA p;
				p.x = 0;
				p.y = 0;
				p.z = 0;
				center_cord3D(iP, t.nvtx, t.pa, p, 100); // ���������� ��������� ������ ��.
				
				if (distance_Point_to_Segment2(p, sp0, sp1) < dopusk) {
					// ����� ����������� �������.

					icounter++;
					xy_arr[icounter].y_function1 = t.potent[iP];
					xy_arr[icounter].y_function2 = tempfiltr[iP];

					switch (idirectional_for_XY_Plot) {
					case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
						xy_arr[icounter].x_argument = p.x;
						break;
					case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: // XZ
						xy_arr[icounter].x_argument = p.y;
						break;
					case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
						xy_arr[icounter].x_argument = p.z;
						break;
					default:
						// YZ
						xy_arr[icounter].x_argument = p.x;
						break;
					}
					
				}
			}
			
			// ���������� �� �����������.
			if (1) {
				std::sort(xy_arr, xy_arr+ icounter+1, compare_XY_PLOT);
			}
			else {
				// BubbleSort
				integer first = 0;
				integer last = icounter;
				integer numberOfPairs = last - first + 1;
				bool swappedElements = true;
				while (swappedElements) {
					numberOfPairs--;
					swappedElements = false;
					for (integer i = first; i <= first + numberOfPairs - 1; i++) {
						if (xy_arr[i].x_argument > xy_arr[i + 1].x_argument) {

							SORT_XY_PLOT temp = xy_arr[i];
							xy_arr[i] = xy_arr[i + 1];
							xy_arr[i + 1] = temp;

							swappedElements = true;
						}
					}
				}
			}

			// ���� ��������� �������� � ����� � ��� �� 
			// ����������, �� �� ��������� ��������� �������,
			// �������� ������ � ����� ������ ��������.
			doublereal pos = xy_arr[0].x_argument;
			integer jposY = 0, jmaximumPos = jposY;
			doublereal maximumTemp = xy_arr[jposY].y_function1;// �����������.
			if (0) {
				for (integer j = 1; j <= icounter; j++) {
					if (pos == xy_arr[j].x_argument) {
						xy_arr[j] = xy_arr[jposY];
					}
					else {
						jposY = j;
						pos = xy_arr[j].x_argument;
					}
				}
			}
			else {
				for (integer j = 1; j <= icounter; j++) {
					if (pos == xy_arr[j].x_argument) {
						if (xy_arr[j].y_function1 > maximumTemp) {
							maximumTemp = xy_arr[j].y_function1;
							jmaximumPos = j;
						}
					}
					else {
						// ������ ����������� �������� �����������
						// � �����������, ��� ���� �������������.
						for (integer j31 = jposY; j31 < j; j31++) {
							xy_arr[j31] = xy_arr[jmaximumPos];
						}
						jposY = j;
						maximumTemp = xy_arr[jposY].y_function1;// �����������.
						jmaximumPos = jposY;
						pos = xy_arr[j].x_argument;
					}
				}
			}

			fprintf(fp, "Xo=%e, Yo=%e, Zo=%e, ", x, y, z);
			switch (idirectional_for_XY_Plot) {
			case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL: // YZ 
				fprintf(fp, "directional=X\n");
				break;
			case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL: // XZ
				fprintf(fp, "directional=Y\n");
				break;
			case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL: // XY
				fprintf(fp, "directional=Z\n");
				break;
			default:
				fprintf(fp, "directional=X\n");
				break;
			}
			fprintf(fp, "position,\ttemperature,\ttemperature_avg\n");
			fprintf(fp1, "VARIABLES= position, temperature\n");
			for (integer j = 0; j <= icounter; j++) {
				fprintf(fp, "%+.16f %+.16f %+.16f\n", xy_arr[j].x_argument, xy_arr[j].y_function1, xy_arr[j].y_function2);
				fprintf(fp1, "%+.16f %+.16f\n", xy_arr[j].x_argument, xy_arr[j].y_function1);
			}
			
			delete[] xy_arr;

			fclose(fp);
			fclose(fp1);
		}
	}
} // xyplot_temp

void animationtecplot360T_3D_part2(integer maxelm, integer ncell,
	FLOW* &f, TEMPER &t, integer flow_interior_count, char* title,
	bool btitle, integer iVar, BLOCK* &b, integer &lb)
{
    // title - ��������� ����.
	// btitle - �������� �� ����� ���������

	FILE *fp=NULL;
    FILE *fp1=NULL; // ����� 1 ��� 3
	
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	
	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	//
#ifdef MINGW_COMPILLER
	int err = 0;
	if (btitle) {
		// �� ������� ���������� ����� � ��������� � ����� ��������:
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
		// �� ��������� ��������� ������������ �����.
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
	errno_t err;

	if (btitle) {
		// �� ������� ���������� ����� � ��������� � ����� ��������:
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
		// �� ��������� ��������� ������������ �����.
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
        
       // int c; // �������� ������
		integer ivarexport=1; // �� ��������� ������ ���� ����������:
		integer i=0; // ������� �����
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
			// ����������� ������ ����� � �������� ����
			// �����������: ������ ���������� �������� ������ ������ � �����:
			if (flow_interior_count>0) {
				// ���� ������ ����. ������ ����� ��������� ���������� ������ ���.
				for (i=0; i<flow_interior_count; i++) if (f[i].bactive) {
					ivarexport=3; // ������� ��� ����������� ������ �������
				}
			}
            
			if (ivarexport==1) {

				// ������ ���������
		        fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

		         // ������ ��� ����������
				switch (iVar) {
	               case TEMP: fprintf(fp, "VARIABLES = x, y, z, Temp\n");  // ���������� ������ ���� ����������
		                       break;
	               case SPEED:  fprintf(fp, "VARIABLES = x, y, z, Speed\n"); // ���������� ������ ��������
		                       break;
	               case PRESS:  fprintf(fp, "VARIABLES = x, y, z, Press\n"); // ���������� ��������
		                       break;
	               case PAM:   fprintf(fp, "VARIABLES = x, y, z, PAM\n"); // ���������� �������� ��������
		                       break;
	            }
		         
#if doubleintprecision == 1
				// ������ ���������� � �����
				fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
				// ������ ���������� � �����
				fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif
				
		
				if (bvery_big_memory) {
					// extended print integer �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (1 && PHYSICAL_MODEL_SWITCH::MESHER_ONLY == steady_or_unsteady_global_determinant) {
						// quolity
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
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
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
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

				    // ������ ���������
		            fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

     				// ������ ��� ����������
				    switch (iVar) {
	                   case TEMP: fprintf(fp, "\nVARIABLES = x, y, z, Temp\n");  // ���������� ������ ���� ����������
		                           break;
	                   case SPEED:  fprintf(fp, "\nVARIABLES = x, y, z, Speed\n"); // ���������� ������ ��������
		                           break;
	                   case PRESS:  fprintf(fp, "\nVARIABLES = x, y, z, Press\n"); // ���������� ��������
		                           break;
	                   case PAM:   fprintf(fp, "\nVARIABLES = x, y, z, PAM\n"); // ���������� �������� ��������
		                           break;
	                }
				}
					
				// ������ ���������� � �����
#if doubleintprecision == 1
				fprintf(fp, "ZONE T=\"%s\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", title, maxelm, ncell);
#else
				fprintf(fp, "ZONE T=\"%s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", title, maxelm, ncell);
#endif
                
				if (bvery_big_memory) {
					// extended print integer �� �������������.

					// ������ x
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
					for (i = 0; i < t.database.maxelm; i++) {
						fprintf(fp, "%+.16f ", t.database.z[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}

					if (1 && PHYSICAL_MODEL_SWITCH::MESHER_ONLY == steady_or_unsteady_global_determinant) {
						// quolity
						for (i = 0; i < t.database.maxelm; i++) {
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
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
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
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
            //fclose(fp1); // �������� �����
            //printf("export tecplot part1 is successfully reading and written...OK.\n");
		}

		// ������ ������ �����
        
        // ������ ���� ���������� ������������ ������.

        if (iVar==TEMP) {
			// ������ �����������
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

		// ������ ����������������� ������� ���� ����������:
		if (ivarexport==3) {

			if (iVar==SPEED) {
				// Speed
                for (i=0;i<maxelm; i++) {
				    if (t.ptr[1][i]>-1) {
					    doublereal svx=f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
					    doublereal svy=f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
					    doublereal svz=f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]*f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
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
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/
			/*
			// VY
            for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/
			/*
			// VZ
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]]);
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
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].diag_coef[VELOCITY_X_COMPONENT][i]);
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
					//fprintf(fp, "%+.16f ", f[t.ptr[1][i]].slau[VELOCITY_X_COMPONENT][i].ap);
				} else fprintf(fp, "%+.16f ", 0.0); 
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/

			/*
			// Mut // ������������ ������������ ��������
			for (i=0;i<maxelm; i++) {
                if (t.ptr[1][i]>-1) {
                   fprintf(fp, "%+.16f ", f[t.ptr[1][i]].potent[MUT][t.ptr[0][i]]);
				} else fprintf(fp, "%+.16f ", 0.0);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");
			*/

			/*
			// ��� ������� ������ ������������� ��������� ����������� �������.
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
				// ����������� ������� ����� � �������� ����
				while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				fclose(fp1); // �������� �����
				//printf("export tecplot part1 is successfully reading and written...OK.\n");
			}
		}
		*/

		//if (bvery_big_memory) 
		{
			// ������ ���������� � ���������� �����
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
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
						// ������������ ������� ���� � ������ ��� � ������.
						//fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
						//fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
						fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


					}

				}
				else {
					//printf("fluid plot\n");
					//getchar();

					if ((ionly_solid_visible == WHAT_VISIBLE_OPTION::ONLY_SOLID_BODY_VISIBLE) && (flow_interior > 0))
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

						integer inode2W = t.neighbors_for_the_internal_node[W_SIDE][0][inode1];
						integer inode3W = t.neighbors_for_the_internal_node[W_SIDE][0][inode4];
						integer inode6W = t.neighbors_for_the_internal_node[W_SIDE][0][inode5];
						integer inode7W = t.neighbors_for_the_internal_node[W_SIDE][0][inode8];


						integer inode5B = t.neighbors_for_the_internal_node[B_SIDE][0][inode1];
						integer inode6B = t.neighbors_for_the_internal_node[B_SIDE][0][inode2];
						integer inode7B = t.neighbors_for_the_internal_node[B_SIDE][0][inode3];
						integer inode8B = t.neighbors_for_the_internal_node[B_SIDE][0][inode4];



						integer inode3S = t.neighbors_for_the_internal_node[S_SIDE][0][inode2];
						integer inode4S = t.neighbors_for_the_internal_node[S_SIDE][0][inode1];
						integer inode7S = t.neighbors_for_the_internal_node[S_SIDE][0][inode6];
						integer inode8S = t.neighbors_for_the_internal_node[S_SIDE][0][inode5];


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

						// ������������ ������ ������� ����.
						// ������������� ==-1 ������� � ��� ��� ���� ����������� ������ ������� ���� � �� ����� ������� �������������.
						if (((t.ptr[1][inode1] == -1) && (t.ptr[1][inode2] == -1) && (t.ptr[1][inode3] == -1) && (t.ptr[1][inode4] == -1) &&
							(t.ptr[1][inode5] == -1) && (t.ptr[1][inode6] == -1) && (t.ptr[1][inode7] == -1) && (t.ptr[1][inode8] == -1)))
						{

							if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

#if doubleintprecision == 1
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
								fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
								fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
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
							// ������������ ������� ���� � ��������.
							// fprintf(fp, "%lld %lld %lld %lld ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%lld %lld %lld %lld\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#else
							// ������������ ������� ���� � ��������.
							// fprintf(fp, "%d %d %d %d ", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i]);
							// fprintf(fp, "%d %d %d %d\n", t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
							fprintf(fp, "%d %d %d %d %d %d %d %d\n", t.database.nvtxcell[0][i], t.database.nvtxcell[1][i], t.database.nvtxcell[2][i], t.database.nvtxcell[3][i], t.database.nvtxcell[4][i], t.database.nvtxcell[5][i], t.database.nvtxcell[6][i], t.database.nvtxcell[7][i]);
#endif


						}

					}
				}
			}
		}


		fclose(fp); // �������� �����
		//printf("export tecplot is successfully written...OK.\n");
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

} // �������� ������ ��������

// �������� ����� ������� � ��������� tecplot 360
// ��-�� ����������� tecplot 360 �� ������ ����� ����� ��������� � ��
// �� ����� ����� ��������� ������� �� ���. ���� ��������� ����������� ��������� �������, ��
// ��� ���� ����� �������� �������� ���������� ������ �� ������ �� ������ �������.
void animationtecplot360T_3D_part2all(integer maxelm, integer ncell,
	FLOW* &f, TEMPER &t, integer flow_interior_count, char* title, bool btitle,
	BLOCK* &b, integer &lb)
{
	// �������� ������� ������� �����������.
	bool bTEMP=true;
	bool bSpeed=true;
	bool bPressure=true;
	bool bPAM=true;

	// ��������������� ������ �������� � ������ �����.
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

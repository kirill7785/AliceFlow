// ���� my_export_tecplot3.c �������� �����������
// ������������� � ��������� tecplot360

#ifndef MY_EXPORT_TECPLOT3_C
#define MY_EXPORT_TECPLOT3_C 1

//#include "windows.h" // ��� ������� WinExec
#include <string.h>

// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
void exporttecplotxy360_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal** potent, doublereal **rhie_chow)
{
	FILE *fp;
	errno_t err;
	// �������� ����� ��� ������.
	if ((err = fopen_s( &fp, "ALICEFLOW0_03.PLT", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

		// ������ ��� ����������
		fprintf(fp, "VARIABLES = x, y, z, Vx, Vy, Vz, Mag, Pressure , Normal , PAm, Zero\n");

		// ������ ���������� � �����
		//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
        //if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
        fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);

		integer i=0; // �������� 
		integer j=0; // ����� for

		// ������ x
	    for (i=0; i<maxelm; i++) {	
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x));
			if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// ������ y
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].y+pa[nvtx[2][i]-1].y));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

        // ������ z
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		
		// ������ �������������� Vx ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[VX][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

        // ������ ������������ Vy ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[VY][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

        // ������ Vz ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[VZ][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");


		// ������ ������ ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", sqrt(potent[VX][i]*potent[VX][i]+potent[VY][i]*potent[VY][i]+potent[VZ][i]*potent[VZ][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// ������ ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[PRESS][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");
        
		// ������ Normal
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", rhie_chow[0][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// ������ PAm
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ",potent[PAM][i] ); // rhie_chow[1][i]
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// ������ Rhi-Chow
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", rhie_chow[2][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// ������ ���������� � ���������� �����
		for (i=0;i<ncell; i++) {
			fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
            fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
		}

		fclose(fp); // �������� �����
		printf("file is successfully written...OK.\n");
       // WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
		  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);
	}
} // exporttecplotxy360_3D

// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
void exporttecplotxy360T_3D(integer maxelm, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, doublereal* potent)
{
	FILE *fp;
	errno_t err;
	// �������� ����� ��� ������.
	if ((err = fopen_s( &fp, "ALICEFLOW0_03.PLT", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");

		// ������ ��� ����������
		fprintf(fp, "VARIABLES = x, y, z, Temp\n");

		// ������ ���������� � �����
        fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);

		integer i=0; // �������� 
		integer j=0; // ����� for

		// ������ x
	    for (i=0; i<maxelm; i++) {	
			fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x));
			if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// ������ y
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].y+pa[nvtx[2][i]-1].y));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

        // ������ z
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		
		// ������ �����������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// ������ ���������� � ���������� �����
		for (i=0;i<ncell; i++) {
			fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
            fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
		}

		fclose(fp); // �������� �����
		printf("file is successfully written...OK.\n");
       // WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
		  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);
	}
} // exporttecplotxy360T_3D

// ���������� ������ ����� ��������� ���������� 19N ����������� ������.
// ����� ����� ����������� ����������� ������������ ��������� �����.

// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 1 � 3.
void exporttecplotxy360T_3D_part1and3(integer maxelm, integer maxbound, bool bextendedprint, integer ncell, 
									  integer** nvtx, integer** nvtxcell, TOCHKA* pa,
									  BOUND* border_neighbor, integer ivarexport)
{
	
	if (bvery_big_memory) {

		database.maxelm = maxelm;
		database.ncell = ncell;

		// extended printeger �� �������������.

		database.x = (doublereal*)malloc(maxelm*sizeof(doublereal));
		database.y = (doublereal*)malloc(maxelm*sizeof(doublereal));
		database.z = (doublereal*)malloc(maxelm*sizeof(doublereal));
		database.nvtxcell = new int*[8];
		for (integer i = 0; i < 8; i++) {
			database.nvtxcell[i] = new integer[ncell];
		}

		// ������ x
		for (integer i = 0; i < maxelm; i++) {
			database.x[i]= 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x);
		}
		// ������ y
		for (integer i = 0; i < maxelm; i++) {
			database.y[i] = 0.5*(pa[nvtx[0][i] - 1].y + pa[nvtx[2][i] - 1].y);
		}

		// ������ z
		for (integer i = 0; i < maxelm; i++) {
			database.z[i]= 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z);
		}

		// ������ ���������� � ���������� �����
		for (integer i = 0; i < ncell; i++) {
			for (integer j = 0; j < 8; j++) {
				database.nvtxcell[j][i] = nvtxcell[j][i];
			}
		}

	}
	else {

		database.maxelm = 0;
		database.ncell = 0;
		//������ � ���� ����� ���������.
		database.x = NULL;
		database.y = NULL;
		database.z = NULL;
		database.nvtxcell = NULL;

		// ����������� ������
		// ��� ����������� ������ �� �������� ����� � ��������� ����.
		// bextendedprint=true; ����������� ������. 

		// ivarexport == 1 ���������� ������ ���� ����������,
		// ivarexport == 2 ���������� ������ �������������,
		// ivarexport == 3 ���������� � ���� ���������� � �������������.

		FILE *fp;
		errno_t err;
		// �������� ����� ��� ������:
		// ���� ������� �� ��� ������: 
		// 1 � 3 ����� ������������ �����
		// ������ ����� � ������������ ������� ������������
		// ����� �������. ����� ���������� ������ ����� ������� � �����
		// ���������� ������ ������������ ����������� ������.
		// �������� ������ 19N.

		// ������ ������ 1 � 3
		if ((err = fopen_s(&fp, "ALICEFLOW0_06_temp_part1.txt", "w")) != 0) {
			printf("Create File temp part1 Error\n");
			getchar();

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

			// ������ ���������� � �����
			fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
			*/

			integer i = 0; // �������� 
			integer j = 0; // ����� for

			// ������ x
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			if (bextendedprint) {
				for (i = 0; i < maxbound; i++) {
					switch (border_neighbor[i].Norm) {// ��������� ���������� ������� � �������
					case E: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI] - 1].x);
						break;
					case W: fprintf(fp, "%+.16f ", pa[nvtx[1][border_neighbor[i].iI] - 1].x);
						break;
					case N: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					case S: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					case T: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
						break;
					case B: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].x + pa[nvtx[1][border_neighbor[i].iI] - 1].x));
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
					case E: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					case W: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					case N: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI] - 1].y);
						break;
					case S: fprintf(fp, "%+.16f ", pa[nvtx[2][border_neighbor[i].iI] - 1].y);
						break;
					case T: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
						break;
					case B: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].y + pa[nvtx[2][border_neighbor[i].iI] - 1].y));
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
					case E: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case W: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case N: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case S: fprintf(fp, "%+.16f ", 0.5*(pa[nvtx[0][border_neighbor[i].iI] - 1].z + pa[nvtx[4][border_neighbor[i].iI] - 1].z));
						break;
					case T: fprintf(fp, "%+.16f ", pa[nvtx[0][border_neighbor[i].iI] - 1].z);
						break;
					case B: fprintf(fp, "%+.16f ", pa[nvtx[4][border_neighbor[i].iI] - 1].z);
						break;
					}
					if (i % 10 == 0) fprintf(fp, "\n");
				}
			}

			fprintf(fp, "\n");

			fclose(fp); // �������� �����
			printf("export tecplot temperature part1 is successfully written...OK.\n");

		}

		if ((err = fopen_s(&fp, "ALICEFLOW0_06_temp_part3.txt", "w")) != 0) {
			printf("Create File temp part3 Error\n");
			getchar();
		}
		else {

			integer i = 0;
			// ������ ���������� � ���������� �����
			for (i = 0; i < ncell; i++) {
				fprintf(fp, "%d %d %d %d ", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
				fprintf(fp, "%d %d %d %d\n", nvtxcell[4][i], nvtxcell[5][i], nvtxcell[6][i], nvtxcell[7][i]);
			}

			fclose(fp); // �������� �����
			printf("export tecplot temperature part3 is successfully written...OK.\n");
		}
	}

} // exporttecplotxy360T_3D_part1and3

// �������� ����������� �����
// ������� ���������� ������� � ��������� tecplot360
// ����� 2.
void exporttecplotxy360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate, bool bextendedprint)
{
    // ianimate - ����� ����������� � ����� ����� ��� ��������.
	bool bprintmessage=false;

	FILE *fp;
    FILE *fp1; // ����� 1 ��� 3
	errno_t err;
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	
	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	// 
	if ((err = fopen_s( &fp, "ALICEFLOW0_07_temp.PLT",  "w")) != 0) {
		printf("Create File temp Error\n");
		getchar();

	}
	else {
        
        char c; // �������� ������
		integer ivarexport=1; // �� ��������� ������ ���� ����������:
		integer i=0; // ������� �����

		bool bOk = true;
		if (!bvery_big_memory) {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
				printf("Open File temp part1 Error\n");
				system("pause");
				bOk = false;

			}
		}
	    if (bOk){
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
		         fprintf(fp, "VARIABLES = x, y, z, Temp, Lam\n");  // ���������� ������ ���� ����������
			
				 // ������ ���������� � �����
				 if (bextendedprint) {
					 fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm+t.maxbound, ncell);
				 }
				 else {
                    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				 }
		
				 if (bvery_big_memory) {
					 // extended printeger �� �������������.

					 // ������ x
					 for (i = 0; i < database.maxelm; i++) {
						 fprintf(fp, "%+.16f ",database.x[i]);
						 if (i % 10 == 0) fprintf(fp, "\n");
					 }
					 // ������ y
					 for (i = 0; i < database.maxelm; i++) {
						 fprintf(fp, "%+.16f ", database.y[i]);
						 if (i % 10 == 0) fprintf(fp, "\n");
					 }
					 // ������ z
					 for (i = 0; i < database.maxelm; i++) {
						 fprintf(fp, "%+.16f ", database.z[i]);
						 if (i % 10 == 0) fprintf(fp, "\n");
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
				if (bextendedprint) {
					 fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z\n");
				}
				else {
                    fprintf(fp, "\nVARIABLES = x, y, z, Temp, Lam, Speed, Pressure, PAM, Vx, Vy, Vz, Rho, Mu, Mut, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, heat_flux_x, heat_flux_y, heat_flux_z\n");
				}
					
				// ������ ���������� � �����
				if (bextendedprint) {
                    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm+t.maxbound, ncell);
				}
				else {
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
				}

				if (bvery_big_memory) {
					// ������ x
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.x[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ y
					for (i = 0; i < database.maxelm; i++) {
						fprintf(fp, "%+.16f ", database.y[i]);
						if (i % 10 == 0) fprintf(fp, "\n");
					}
					// ������ z
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

		// ������ ����������������� ������� ���� ����������:
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
				// �������� � ������ ���� ��� ������������ �������� ������ ����, ���� ������ ���� ���������� �������.
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

			// Mut // ������������ ������������ ��������
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

			// ��� ������� ������ ������������� ��������� ����������� �������.
			// ��� ������������� ���������� �� ������ Distance_Wall.
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

			

			// Curl // ������������ - ������ ������ ��������.
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

			// ������� ����������� �� ��������� �������� !!!.

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
					t.border_neighbor,  Tx, Ty, Tz);
			}

			for (i=0; i<t.maxelm; i++) {
				// ������ ��������� ����.
				green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
					t.neighbors_for_the_internal_node, t.maxelm, true, 
					t.border_neighbor,Tx, Ty, Tz);
			}
			

			// ���������� � ����.

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

			// ������������ ����������� ������.
			delete Tx;
			delete Ty;
			delete Tz;

		}

		if (bvery_big_memory) {
			// ������ ���������� � ���������� �����
			for (i = 0; i < database.ncell; i++) {
				fprintf(fp, "%d %d %d %d ", database.nvtxcell[0][i], database.nvtxcell[1][i], database.nvtxcell[2][i], database.nvtxcell[3][i]);
				fprintf(fp, "%d %d %d %d\n", database.nvtxcell[4][i], database.nvtxcell[5][i], database.nvtxcell[6][i], database.nvtxcell[7][i]);
			}
		}
		else {
			if ((err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
				printf("Open File temp part3 Error\n");
				getchar();

			}
			else {
				// ����������� ������� ����� � �������� ����
				while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				fclose(fp1); // �������� �����
				if (bprintmessage) {
					printf("export tecplot part1 is successfully reading and written...OK.\n");
				}
			}
		}

		fclose(fp); // �������� �����
		if (bprintmessage) {
			printf("export tecplot is successfully written...OK.\n");
		}
		else printf("export tecplot 360... "); // �������� ��������� ��� �������� �� ����� ������.
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
		doublereal x=0.0e-3, y=0.0e-3, z=0.0e-3; // ����� ����� ������� �������� �����
		integer ifi=0, iPf=0;
		integer iplane=XZ; // ��������� ���������������� �����.
		for (integer i=0; i<flow_interior; i++) {
			for (integer iP=0; iP<fglobal[i].maxelm; iP++) {
				center_cord3D(iP, fglobal[ifi].nvtx, fglobal[ifi].pa, p); // ���������� ��������� ������ ��.
				dist=sqrt(fabs(x-p.x)*fabs(x-p.x)+fabs(y-p.y)*fabs(y-p.y)+fabs(z-p.z)*fabs(z-p.z));
				if (dist<epsilon) {
					epsilon=dist;
					ifi=i;
					iPf=iP;
				}
			}
		}
		// ��������� � ������ 
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
		doublereal dx=0.0, dy=0.0, dz=0.0;// ����� �������� ������������ ������
	    volume3D(iPf, fglobal[ifi].nvtx, fglobal[ifi].pa, dx, dy, dz);
        center_cord3D(iPf, fglobal[ifi].nvtx,fglobal[ifi].pa, p); // ���������� ��������� ������ ��.
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
					// ����� ������������������ ����� ��� �������.
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

// ���������� ������� ����������� ����� ����� 
void xyplot_temp(TEMPER &t, doublereal* tempfiltr) {

	// tempfiltr - ������������ ���������� ������������� �����������.
	// ��� ������������ ����� xyplotT.txt.

	FILE *fp;
	errno_t err;
	if ((err = fopen_s( &fp, "xyplotT.txt",  "w")) != 0) {
		printf("Create File xyplot Error\n");
		getchar();

	}
	else {

		FILE* fp1;
		errno_t err1;
		if ((err1 = fopen_s(&fp1, "xyplotT1.PLT", "w")) != 0) {
			printf("Create File xyplot Error\n");
			getchar();

		}
		else {

			// �������� ! ��������� ������� ����� ����� ������� ����� ��������� ����� � ��������� ������� ������ ����� ����� ���������������.
			TOCHKA p;
			doublereal epsilon = 1e30;
			doublereal dist;
			doublereal x = -13.25e-6, y = 0.0e-6, z = 0.0e-6; // ����� ����� ������� �������� �����
			integer iPf = 0;
			integer iplane = XY; // ��������� ���������������� �����.

			for (integer iP = 0; iP < t.maxelm; iP++) {
				center_cord3D(iP, t.nvtx, t.pa, p); // ���������� ��������� ������ ��.
				dist = sqrt(fabs(x - p.x) * fabs(x - p.x) + fabs(y - p.y) * fabs(y - p.y) + fabs(z - p.z) * fabs(z - p.z));
				if (dist < epsilon) {
					epsilon = dist;
					iPf = iP;
				}
			}

			// ��������� � ������ 
			switch (iplane) {
			case XY: while (t.neighbors_for_the_internal_node[B][iPf].iNODE1 < t.maxelm) iPf = t.neighbors_for_the_internal_node[B][iPf].iNODE1; break;
			case XZ: while (t.neighbors_for_the_internal_node[S][iPf].iNODE1 < t.maxelm) iPf = t.neighbors_for_the_internal_node[S][iPf].iNODE1; break;
			case YZ: while (t.neighbors_for_the_internal_node[W][iPf].iNODE1 < t.maxelm) iPf = t.neighbors_for_the_internal_node[W][iPf].iNODE1; break;
			}

			integer G;
			switch (iplane) {
			case XY: G = T;  break;
			case XZ: G = N;  break;
			case YZ: G = E;  break;
			}


			fprintf(fp, "position,\ttemperature,\ttemperature_avg\n");
			fprintf(fp1, "VARIABLE= position temperature\n");
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// ����� �������� ������������ ������
			volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
			center_cord3D(iPf, t.nvtx, t.pa, p); // ���������� ��������� ������ ��.
			switch (iplane) {
			case XY: fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.z - 0.5 * dz, t.potent[t.neighbors_for_the_internal_node[B][iPf].iNODE1],
				tempfiltr[t.neighbors_for_the_internal_node[B][iPf].iNODE1]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.z - 0.5 * dz, t.potent[t.neighbors_for_the_internal_node[B][iPf].iNODE1]);
				break;
			case XZ:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.y - 0.5 * dy, t.potent[t.neighbors_for_the_internal_node[S][iPf].iNODE1],
				tempfiltr[t.neighbors_for_the_internal_node[S][iPf].iNODE1]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.y - 0.5 * dy, t.potent[t.neighbors_for_the_internal_node[S][iPf].iNODE1]);
				break;
			case YZ:  fprintf(fp, "%+.16f %+.16f %+.16f\n",
				p.x - 0.5 * dx, t.potent[t.neighbors_for_the_internal_node[W][iPf].iNODE1],
				tempfiltr[t.neighbors_for_the_internal_node[W][iPf].iNODE1]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.x - 0.5 * dx, t.potent[t.neighbors_for_the_internal_node[W][iPf].iNODE1]);
				break;
			}
			switch (iplane) {
			case XY: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.z, t.potent[iPf],
					tempfiltr[iPf]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.z, t.potent[iPf]);
				if (t.neighbors_for_the_internal_node[T][iPf].iNODE1 >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.z + 0.5 * dz, t.potent[t.neighbors_for_the_internal_node[T][iPf].iNODE1],
						tempfiltr[t.neighbors_for_the_internal_node[T][iPf].iNODE1]);
					fprintf(fp1, "%+.16f %+.16f\n",
						p.z + 0.5 * dz, t.potent[t.neighbors_for_the_internal_node[T][iPf].iNODE1]);

				}
				iPf = t.neighbors_for_the_internal_node[T][iPf].iNODE1;
			} break;
			case XZ: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.y, t.potent[iPf],
					tempfiltr[iPf]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.y, t.potent[iPf]);
				if (t.neighbors_for_the_internal_node[N][iPf].iNODE1 >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.y + 0.5 * dy, t.potent[t.neighbors_for_the_internal_node[N][iPf].iNODE1],
						tempfiltr[t.neighbors_for_the_internal_node[N][iPf].iNODE1]);
					fprintf(fp1, "%+.16f %+.16f\n",
						p.y + 0.5 * dy, t.potent[t.neighbors_for_the_internal_node[N][iPf].iNODE1]);
				}
				/*
				// ����� ������������������ ����� ��� �������.
				printf("iPf=%d\n",iPf);
				if (fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1>=fglobal[ifi].maxelm) {
					printf("iPffinish=%d\n",fglobal[ifi].neighbors_for_the_internal_node[N][iPf].iNODE1);
					getchar();
				}
				*/
				iPf = t.neighbors_for_the_internal_node[N][iPf].iNODE1;
			} break;
			case YZ: while (iPf < t.maxelm) {
				center_cord3D(iPf, t.nvtx, t.pa, p);
				fprintf(fp, "%+.16f %+.16f %+.16f\n",
					p.y, t.potent[iPf],
					tempfiltr[iPf]);
				fprintf(fp1, "%+.16f %+.16f\n",
					p.y, t.potent[iPf]);
				if (t.neighbors_for_the_internal_node[E][iPf].iNODE1 >= t.maxelm) {
					volume3D(iPf, t.nvtx, t.pa, dx, dy, dz);
					fprintf(fp, "%+.16f %+.16f %+.16f\n",
						p.x + 0.5 * dx, t.potent[t.neighbors_for_the_internal_node[E][iPf].iNODE1],
						tempfiltr[t.neighbors_for_the_internal_node[E][iPf].iNODE1]);
					fprintf(fp1, "%+.16f %+.16f\n",
						p.x + 0.5 * dx, t.potent[t.neighbors_for_the_internal_node[E][iPf].iNODE1]);
				}
				iPf = t.neighbors_for_the_internal_node[E][iPf].iNODE1;
			} break;
			}
		}
		fclose(fp);
		fclose(fp1);
	}
} // xyplot_temp

void animationtecplot360T_3D_part2(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, char* title, bool btitle, integer iVar)
{
    // title - ��������� ����.
	// btitle - �������� �� ����� ���������

	FILE *fp=NULL;
    FILE *fp1; // ����� 1 ��� 3
	errno_t err;
	// �������� ����� ��� ������:
	// ���� ������� �� ��� ������: 
	// 1 � 3 ����� ������������ �����
	// ������ ����� � ������������ ������� ������������
	// ����� �������. ����� ���������� ������ ����� ������� � �����
	// ���������� ������ ������������ ����������� ������.
	// �������� ������ 19N.

	
	// ������ ������ 1 � 3 � ������ ���� ��� ������ � �������� ����.
	//
	if (btitle) {
		// �� ������� ���������� ����� � ��������� � ����� ��������:
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
		// �� ��������� ��������� ������������ �����.
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
        
        char c; // �������� ������
		integer ivarexport=1; // �� ��������� ������ ���� ����������:
		integer i=0; // ������� �����

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part1.txt", "r")) != 0) {
		    printf("Open File temp part1 Error\n");
		    getchar();

	    }
	    else {
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
		         
			
				 // ������ ���������� � �����
                fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
		

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

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
                fprintf(fp, "ZONE T=\"%s\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", title,  maxelm, ncell);

				while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);

						
			}
            fclose(fp1); // �������� �����
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

        if ((err = fopen_s( &fp1, "ALICEFLOW0_06_temp_part3.txt", "r")) != 0) {
		    printf("Open File temp part3 Error\n");
		    getchar();

	    }
	    else {
			// ����������� ������� ����� � �������� ����
			while ( (c=fgetc(fp1))!=EOF) fputc(c,fp);
            fclose(fp1); // �������� �����
            //printf("export tecplot part1 is successfully reading and written...OK.\n");
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
void animationtecplot360T_3D_part2all(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, char* title, bool btitle)
{
	// �������� ������� ������� �����������.
	bool bTEMP=true;
	bool bSpeed=true;
	bool bPressure=true;
	bool bPAM=true;

	// ��������������� ������ �������� � ������ �����.
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

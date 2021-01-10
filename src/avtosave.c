// ���� avtosave.c �������������� � ������������� �������
// � ����������� �����.
// ������ 16 ������� 2011 �.

#ifndef AVTO_SAVE_C
#define AVTO_SAVE_C 1

// �������������� ������������ ������� ������� � ���� avtosave.txt
// ��� ������������� ����� ����� ����������.
void avtosave(FLOW* &f, TEMPER &t, integer flow_interior, integer* &inumber_iteration_SIMPLE, doublereal* &continity_start) {

	// flow_interior - ����� ����� ����������� FLUID ���.
	// inumber_iteration_SIMPLE - ����� �������� SIMPLE ��������� �� ������ ������� ���������� ������ � ���� (��� ������ ������ ����).
	// continity_start - ��������� �������� ������� continity (��� ������ ������ ����).


	FILE *fp_avtosave=NULL; // ���� � ������� ����� ������������ ����������� ����
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp_avtosave=fopen64("avtosave.txt", "w");
	if (fp_avtosave == NULL) {
		err = 1;
	}
#else
	errno_t err = 0;
	err = fopen_s(&fp_avtosave, "avtosave.txt", "w");
#endif
	
	if (((err) !=0)||(fp_avtosave == NULL)) {
		// ���� ����������� ��� ������ ����������,
		// ���� ������ ����� ������ �� ���� �� �� ����� ������.
	   printf("Create File avtosave.txt Error\n");
       //getchar();
	   system("pause");
       exit(0);
	}
	else {

#if doubleintprecision == 1
		fprintf(fp_avtosave, "%lld %lld ", t.maxelm + t.maxbound, flow_interior);
		integer i = 0; // ������� ������ FLUID ���
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%lld ", f[i].maxelm + f[i].maxbound);
		fprintf(fp_avtosave, "\n");

		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%lld ", inumber_iteration_SIMPLE[i]);
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%+.16f ", continity_start[i]);
		fprintf(fp_avtosave, "\n"); // ������� � ����� ������
#else
		fprintf(fp_avtosave, "%d %d ", t.maxelm + t.maxbound, flow_interior);
		integer i = 0; // ������� ������ FLUID ���
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%d ", f[i].maxelm + f[i].maxbound);
		fprintf(fp_avtosave, "\n");

		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%d ", inumber_iteration_SIMPLE[i]);
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%+.16f ", continity_start[i]);
		fprintf(fp_avtosave, "\n"); // ������� � ����� ������
#endif

		

		// TEMP
		for (i=0; i<t.maxelm+t.maxbound; i++) {
           fprintf(fp_avtosave,"%+.16f ",t.potent[i]);
		   if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		}
		if ((t.maxelm+t.maxbound)%20 != 0) fprintf(fp_avtosave,"\n");

		integer j=0;
		for (j=0; j<flow_interior; j++) {
			// VX
			for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VELOCITY_X_COMPONENT][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// VY
            for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VELOCITY_Y_COMPONENT][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// VZ
			for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VELOCITY_Z_COMPONENT][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// PRESSURE
			for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[PRESS][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");

		} // ����� ����� �� ���� ������ �����

		fclose(fp_avtosave); // �������� ����� ��� ������ ������������ ������� �������.
	}


} // avtosave

// ������������� ����� ����� ����������:
// ���������� �������� ������� ������� �� ����� avtosave.txt 
// ��� ��� ������� � ������������.
void avtoreadvalue(FLOW* &f, TEMPER &t, integer flow_interior,
	integer* &inumber_iteration_SIMPLE, doublereal* &continity_start, bool &breadOk,
	BLOCK* &b, integer &lb, SOURCE* &s_loc, integer &ls, WALL* &w, integer &lw) {

	// flow_interior - ����� ����� ����������� FLUID ���.
	// inumber_iteration_SIMPLE - ����� �������� SIMPLE ��������� �� ������ ������� ���������� ������ � ���� (��� ������ ������ ����).
	// continity_start - ��������� �������� ������� continity (��� ������ ������ ����).

	FILE *fp_avtosave=NULL; // ���� � ������� �������� ���������� ������� ��������.
	
#ifdef MINGW_COMPILLER
	int err1 = 0;
	fp_avtosave=fopen64("avtosave.txt", "r");
	if (fp_avtosave == NULL) err1 = 1;
#else
	errno_t err1 = 0;
	err1 = fopen_s(&fp_avtosave, "avtosave.txt", "r");
#endif

	if ((err1)!=0) {
		// ����� �������������� ���������� ������ ���,
		// ������ ���� ������� �� �����.
		printf("No input File avtosave.txt \n");
		breadOk=false;
		//getchar();
	}
	else if (fp_avtosave!=NULL)
	{
		bool bcontinue=true; // ���������� �� ���������� �����

		//float fin=0.0;
		doublereal fin=0.0;
		integer din=0;

		// ������� ��������� �������� ������������ ������:
		// premeshin.txt � avtosave.txt
		// bcontinue==true ���� �������� ������ �������.
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
		fscanf(fp_avtosave, "%lld", &din);
#else
		fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
		fscanf_s(fp_avtosave, "%lld", &din);
#else
		fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif
        
		if (din != (t.maxelm+t.maxbound)) bcontinue=false; // �������� �� ���������� ����� �� � ������ ����������������

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
		fscanf(fp_avtosave, "%lld", &din);
#else
		fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
		fscanf_s(fp_avtosave, "%lld", &din);
#else
		fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif

		
		if (din != flow_interior) bcontinue=false; // �������� �� ���������� ������ FLUID ���

		integer i=0; // ������� FLUID ���
		for (i=0; i<flow_interior; i++) {

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
			fscanf(fp_avtosave, "%lld", &din);
#else
			fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
			fscanf_s(fp_avtosave, "%lld", &din);
#else
			fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif

			
			if (din != (f[i].maxelm+f[i].maxbound)) bcontinue=false; // �������� �� ���������� �� � FLUID ���� ����� i.
		}

		if (!bcontinue) {
			// ���� premeshin.txt �� ������������� ����� avtosave.txt
			printf("premeshin.txt does not correspond to avtosave.txt...\n");
			// ���������� ������� ������
			printf("computation begins a new! Please press any key to continue...\n");
			breadOk=false;
			//getchar();
			system("pause");
		}
		else {

			for (i=0; i<flow_interior; i++) {

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
				fscanf(fp_avtosave, "%lld", &din);
#else
				fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
				fscanf_s(fp_avtosave, "%lld", &din);
#else
				fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif

				
				inumber_iteration_SIMPLE[i]=din;
			}
			for (i=0; i<flow_interior; i++) {

#ifdef MINGW_COMPILLER
				fscanf(fp_avtosave, "%lf", &fin);
#else
				fscanf_s(fp_avtosave, "%lf", &fin);
#endif

				
				continity_start[i]=fin;
			}

			integer j=0; // ������� ��.
			// ������ ���� ����������
			printf("read TEMP...\n");
			for (j=0; j<(t.maxelm+t.maxbound); j++) {

#ifdef MINGW_COMPILLER
				fscanf(fp_avtosave, "%lf", &fin);
#else
				fscanf_s(fp_avtosave, "%lf", &fin);
#endif

               
				t.potent[j]=fin;
			}
			// ���� �� ���� ������ �����:
			for (i=0; i<flow_interior; i++) {
#if doubleintprecision == 1
				printf("read FLUID%lld...\n", i);
#else
				printf("read FLUID%d...\n", i);
#endif
				
				// VX
				printf("read VX...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VELOCITY_X_COMPONENT][j]=fin;
					f[i].potent[VXCOR][j]=fin;
				}
				// VY
				printf("read VY...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VELOCITY_Y_COMPONENT][j]=fin;
					f[i].potent[VYCOR][j]=fin;
				}
				// VZ
				printf("read VZ...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VELOCITY_Z_COMPONENT][j]=fin;
					f[i].potent[VZCOR][j]=fin;
				}
				// PRESSURE
				printf("read PRESSURE...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[PRESS][j]=fin;
				}

				// ���������� ���������� ��������:
				for (j=0; j<f[i].maxelm; j++) {
					// ��������� ��������� ��� ���������� ��.
				    green_gauss(j, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, false,f[i], f[i].border_neighbor, t.ilevel_alice);
				}
				for (j=0; j<f[i].maxelm; j++) {
					// ��������� ��������� ��� ��������� ��.
				    green_gauss(j, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, true, f[i], f[i].border_neighbor, t.ilevel_alice);
				}

				// ���������� ���������� ��������:
	            // �� ������ ������������������ ���� ��������.
            	// ��������� �������� ����������� ��� ���������� �������� ���-��� 1983�.
				for (integer j_1 = 0; j_1 < f[i].maxelm; j_1++) {
					// ��������� �������� ��� ���������� ��.
					green_gaussPRESS(j_1, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, false, f[i].border_neighbor, ls, lw, w, f[i].bLR1free, t.ilevel_alice, f[i].ptr, f[i].volume);
				}
				for (integer j_1 = 0; j_1 < f[i].maxelm; j_1++) {
					// ��������� �������� ��� ��������� ��.
					green_gaussPRESS(j_1, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, true, f[i].border_neighbor, ls, lw, w, f[i].bLR1free, t.ilevel_alice, f[i].ptr, f[i].volume);
				}



				// �� ������ ������ ������������ ������������ �������� ����� ����.
	            // ���������� S ���������� ������� ���������-���������� ��� ����
	            // ���������� � ��������� ����������� ������� ��������� �������.
	            #pragma omp parallel for shared (f) schedule (guided)
		        for (integer j_1=0; j_1<(f[i].maxelm+f[i].maxbound); j_1++) {
			        // �� ������ ������������ ������� ��. user_manual.
			        doublereal sum=0.0;
			        sum+=2.0*f[i].potent[GRADXVX][j_1]*f[i].potent[GRADXVX][j_1];
			        sum+=2.0*f[i].potent[GRADYVY][j_1]*f[i].potent[GRADYVY][j_1];
			        sum+=2.0*f[i].potent[GRADZVZ][j_1]*f[i].potent[GRADZVZ][j_1];
			        sum+=(f[i].potent[GRADYVX][j_1]+f[i].potent[GRADXVY][j_1])*(f[i].potent[GRADYVX][j_1]+f[i].potent[GRADXVY][j_1]);
			        sum+=(f[i].potent[GRADZVX][j_1]+f[i].potent[GRADXVZ][j_1])*(f[i].potent[GRADZVX][j_1]+f[i].potent[GRADXVZ][j_1]);
			        sum+=(f[i].potent[GRADYVZ][j_1]+f[i].potent[GRADZVY][j_1])*(f[i].potent[GRADYVZ][j_1]+f[i].potent[GRADZVY][j_1]);
					// ��������� ��������� ����� ������� ����������� ��������� �������������.
		            // �������� ��� ����� �������� �����������.
					sum-=(2.0/3.0)*(f[i].potent[GRADXVX][j_1]+f[i].potent[GRADYVY][j_1]+f[i].potent[GRADZVZ][j_1])*(f[i].potent[GRADXVX][j_1]+f[i].potent[GRADYVY][j_1]+f[i].potent[GRADZVZ][j_1]); // ������� ��������� � ��������������/������������
			        f[i].SInvariantStrainRateTensor[j_1]=sqrt(fmax(0.0,sum));

                    // ����� (������ ������ ��������).
					// �� ����� ��� ���� ����� ��� � ������ ������� ��������.
		            sum=0.0;
		            sum+=(f[i].potent[GRADYVZ][j_1]-f[i].potent[GRADZVY][j_1])*(f[i].potent[GRADYVZ][j_1]-f[i].potent[GRADZVY][j_1]); // ��������� !
		            sum+=(f[i].potent[GRADZVX][j_1]-f[i].potent[GRADXVZ][j_1])*(f[i].potent[GRADZVX][j_1]-f[i].potent[GRADXVZ][j_1]);
		            sum+=(f[i].potent[GRADXVY][j_1]-f[i].potent[GRADYVX][j_1])*(f[i].potent[GRADXVY][j_1]-f[i].potent[GRADYVX][j_1]);
		            f[i].potent[CURL][j_1]=sqrt(sum);
		        }

				breadOk=true; // ���������� ������ �������.

			} // ����� ����� �� ������ �����
			
		    

		}

		

		fclose(fp_avtosave); // �������� ����� �� �������� ����������� ���������� ������� ��������
	}

} // avtoreadvalue


#endif
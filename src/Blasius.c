// ���� Blasius.c 
// 1. ���������� ������� ������������������ ������������ ���� ��� ������ �������� 1908 ���� (Boundary layer thickness);
// 2. ���������� ������������� ���������� ������������ ������ �� �������� (Wall skin friction distribution);
// 3. ���������� ������� ���������� �� ����������� �������� (displacement thickness);
// 4. ���������� ������� �������������� ������������ ���� (Thermal boundary layer distribution);
// 5. ���������� ������������� ����� ��������� �� ����������� ��������.
// begin: 31 ������ 2012 ����.
// end first: 5 ������� 2012 ����.

// ��� ������ �������� ��� ����� ������� �������� ������� 2-3 ����� ���������,
// � ����� 3-5 ���� ���������.

#pragma once
#ifndef MY_BLASIUS_C
#define MY_BLASIUS_C 1

// �������������� ��� ������ �������� 1908 ����.
void boundarylayer_info(FLOW* &f, TEMPER &t, integer flow_interior_count, WALL* w, integer lw) {

	// ����������� ����� ���������� ��� ��������� ������� �������� ����� Re_��=1e+5;

	// �������������� ��� ������ ������� � ��������� � ��������� YZ, 
	// � �������� VX ����� ����.
	// �.�. ��������� ��������� ���������� ��������� �� ��� ���������� ������ ���������
	// ����������� �� ������ ��� X: avgX=0.5*(minX+maxX);

	// �������� ��������� � ���� Y ��� Z==0.0. �������� ���������� ��� Y==0.0 � ������ 1 ����.

	const doublereal length_plate=1.0; // m
	doublereal U_inf=0.0; // �������� ����������� ������ �� �������������.
	doublereal Twall=-1e3; // ����������� �������� ��������
	doublereal Tinf=1e5; // ����������� ������� �� �����.

	// ����������� �������� ����������� ������ �� �������������:
	for (integer i=0; i<lw; i++) {
		if ((!w[i].bsymmetry)&&(!w[i].bpressure)) {
			// ������ �� �������� �� �������� ��������� �� �������� ��������.
			U_inf=fmax(U_inf,w[i].Vy); // �.�. U_inf ��������� �� ����������� � ������������� ������������ ��� OY.
			Twall=fmax(Twall,w[i].Tamb); // ������� ����������� �������� ��������.
			Tinf=fmin(Tinf,w[i].Tamb); // ������� ����������� ������������ �������� ����������� ������.
		}
	}

	// �������� �������������� �������� ������������ ����������� ����������� �������, �.�. Twall>Tinf.

    // ������ � U_inf - �������� �������� ����������� ������ �� �������������.
	doublereal rho_avg=0.0; // ������� ��������� �� ������.
	doublereal mu_avg=0.0; // ������� �� ������ ������������ ��������.
	doublereal lam_avg=0.0; // ������� �� ������ ����������������
	doublereal cp_avg=0.0; // ������� �� ������ �����������
	doublereal volume_default_interior=0.0; // ����� ��������� �������. 
	doublereal VyMAX=-1e30; // �������� �� ��������� ������� ���������� �������� ������������ ��������

	for (integer iP=0; iP<f[0].maxelm; iP++) {
		// ���������� �������� �������� ������������ ������:
	    doublereal dx=0.0, dy=0.0, dz=0.0;// ����� �������� ������������ ������
	    volume3D(iP, f[0].nvtx, f[0].pa, dx, dy, dz);

		rho_avg+=dx*dy*dz*f[0].prop[RHO][iP];
		mu_avg+=dx*dy*dz*f[0].prop[MU_DYNAMIC_VISCOSITY][iP];
		lam_avg+=dx*dy*dz*t.prop[LAM][f[0].ptr[iP]];
		cp_avg+=dx*dy*dz*t.prop[HEAT_CAPACITY][f[0].ptr[iP]];
		volume_default_interior+=dx*dy*dz;

		VyMAX=fmax(VyMAX,f[0].potent[VELOCITY_Y_COMPONENT][iP]);
	}

	rho_avg=rho_avg/volume_default_interior; // ������� �� ������ ��������� �����.
	mu_avg=mu_avg/volume_default_interior; // ������� �� ������ ������������ ��������.
	lam_avg=lam_avg/volume_default_interior; // ������� �� ������ ����������������.
	cp_avg=cp_avg/volume_default_interior; // ������� �� ������ �����������.

	doublereal avg_Re_number=U_inf*length_plate*rho_avg/mu_avg; // ����� ���������� �� ����� ���������.
	doublereal avg_Pr_number=mu_avg*cp_avg/lam_avg; // ����� ��������.
	doublereal avg_Pe_number=avg_Re_number*avg_Pr_number; // ����� �����.

	// ���������� ������������ ������ ����������� �������� � ����� avgX=0.5*(minX+maxX);
	doublereal dx=0.0, dy=0.0, dz=0.0;// ����� �������� ������������ ������
	TOCHKA p; // ���������� ������ �������� ������������ ������.
	integer iP=0;
	doublereal minX, maxX, avgX;
	while (f[0].neighbors_for_the_internal_node[W_SIDE][0][iP]<f[0].maxelm) iP=f[0].neighbors_for_the_internal_node[W_SIDE][0][iP];
    // ���������� �������� �������� ������������ ������:
	volume3D(iP, f[0].nvtx, f[0].pa, dx, dy, dz);
	center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // ���������� ��������� ������ ��.
	minX=p.x-0.5*dx;
	while (f[0].neighbors_for_the_internal_node[E_SIDE][0][iP]<f[0].maxelm) iP=f[0].neighbors_for_the_internal_node[E_SIDE][0][iP];
	volume3D(iP, f[0].nvtx, f[0].pa, dx, dy, dz);
	center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // ���������� ��������� ������ ��.
	maxX=p.x+0.5*dx;
	avgX=0.5*(minX+maxX); // ���������� ����������� ��������� � ������� � ����� ����������� ���������.

	center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // ���������� ��������� ������ ��.
	doublereal mindist=fabs(p.x-avgX);
	integer iPC=iP;
	while (f[0].neighbors_for_the_internal_node[W_SIDE][0][iP]<f[0].maxelm) {
		iP=f[0].neighbors_for_the_internal_node[W_SIDE][0][iP];
		center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // ���������� ��������� ������ ��.
		if (fabs(p.x-avgX)<mindist) {
			iPC=iP;
			mindist=fabs(p.x-avgX);
		}
	}

	// ����������� ����� iPC �������� ��������� � ����������� ��������� � ������� ����� ������������ ������.
	iP=iPC;

	// ���������� ���������� ����������� ������� ������������� �� ����� ��������:
	integer iclength=0;
	while (f[0].neighbors_for_the_internal_node[S_SIDE][0][iP]<f[0].maxelm) iP=f[0].neighbors_for_the_internal_node[S_SIDE][0][iP];

	while (iP<f[0].maxelm) {
		center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // ���������� ��������� ������ ��.
		if ((p.y>0.0) && (p.y<1.0)) iclength++; // �������� ����������� ����� 0.0 m � 1.0 m.
		iP=f[0].neighbors_for_the_internal_node[N_SIDE][0][iP];
	}

	// ����������������� ����������� ����: (Boundary layer thickness).
	doublereal* delta=new doublereal[iclength];
	// ��������������� ���������� �� �������� ������ ��������:
	doublereal* yposition=new doublereal[iclength];
	// ������������ �������� ��������� ����������� ������ �� ����������� ��������:
	// Cx[i]=SInvariantStrainRateTensor[i]/(0.5*rho_avg*U_inf*U_inf); ��. �����.
	// Wall skin friction distribution
	doublereal* Cx=new doublereal[iclength];
	// ������� ����������� ������:
	doublereal avg_Cx=0.0; // integral(Cx[i]*dlength_plate)/length_plate;
	// ������� ���������� �� ����������� ��������. (displacement thickness) 
	// ��. �.�.������, �.�.������ ������������� ������. ������������� ��� VI. ���. 228.
	doublereal* displacement_thickness=new doublereal[iclength];
	// ������������� ����������� ���� (Thermal boundary layer distribution):
	doublereal* deltaT=new doublereal[iclength];
	// ��������� ����� ���������� ��������� �� ������:
	// ����� ��������� ���������� ��������� ������������� ����������� �� ���� ��������� �
	// � ������������� ����������� �� ���� ����������������.
	doublereal* local_Nusselt_number=new doublereal[iclength];


	// ���������� 95% ������� ������������������ ������������ ���� �� ��������:
	doublereal deltascal=0.95; // 95% ����������� ����.
	doublereal deltaTscal=0.05; // 5% �������� ����������� ���� �� �������� ��������.
	iP=iPC;
	while (f[0].neighbors_for_the_internal_node[S_SIDE][0][iP]<f[0].maxelm) iP=f[0].neighbors_for_the_internal_node[S_SIDE][0][iP]; // ��������� � ������.

	integer ilengthcounter=0;
	for (ilengthcounter=0; ilengthcounter<iclength; ilengthcounter++) {
		// �������������.
		delta[ilengthcounter]=0.0; // ����������������� ����������� ���� 
		yposition[ilengthcounter]=0.0; // �������� ��� ��������
		Cx[ilengthcounter]=0.0; // ������������ ��������� ����������� ������ �� ������
        displacement_thickness[ilengthcounter]=0.0; // ������� ����������
		deltaT[ilengthcounter]=0.0; // ������������� ����������� ����
		local_Nusselt_number[ilengthcounter]=0.0; // ��������� ����� ���������.
	}
	ilengthcounter=0;

	while (iP<f[0].maxelm) {
		center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // ���������� ��������� ������ ��.
		if ((p.y>0.0) && (p.y<length_plate)) {
			// �������� ����������� ����� 0.0 m � length_plate m.
			yposition[ilengthcounter]=p.y; // �������� ��� ������� ������������ ����.
			
			integer iBT=iP;
			while (f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]<f[0].maxelm) iBT=f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]; // ��������� � ������ ��������.
			// ���������� ������������� ������������� ���������� ������������ ������ �� ������:
			// ��. �.�.�����  ��������������� ���������� ������., "�������", 1978�. ���. 184.
			Cx[ilengthcounter]=mu_avg*f[0].SInvariantStrainRateTensor[f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]]/(0.5*rho_avg*U_inf*U_inf);

			// ���������� ������� ������������������ ������������ ����.
			doublereal VYB, VYT, zB, zT;
			VYB=f[0].potent[VELOCITY_Y_COMPONENT][f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]];
			VYT=f[0].potent[VELOCITY_Y_COMPONENT][iBT];
			center_cord3D(iBT, f[0].nvtx, f[0].pa, p, B_SIDE); // ���������� ��������� ������ ��.
			volume3D(iBT, f[0].nvtx, f[0].pa, dx, dy, dz);
			zB=p.z-0.5*dz;
			zT=p.z;
			while ((f[0].neighbors_for_the_internal_node[T_SIDE][0][iBT]<f[0].maxelm) && (VYT<deltascal*U_inf)) {
				iBT=f[0].neighbors_for_the_internal_node[T_SIDE][0][iBT]; // ��������� �� �������� ��������������� � ���������.
				VYB=f[0].potent[VELOCITY_Y_COMPONENT][f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]];
			    VYT=f[0].potent[VELOCITY_Y_COMPONENT][iBT];
				zB=zT;
				center_cord3D(iBT, f[0].nvtx, f[0].pa, p, T_SIDE); // ���������� ��������� ������ ��.
				zT=p.z;
			}
			doublereal a, b;
			a=(VYT-VYB)/(zT-zB);
			b=(zT*VYB-VYT*zB)/(zT-zB);
			delta[ilengthcounter]=(deltascal*U_inf-b)/a; // ������� ������� ������������������ ������������ ����.

			while (f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]<f[0].maxelm) iBT=f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]; // ��������� � ������ ��������.
			// ���������� ������� ���������� �� �������� Uoperating (displacement thickness):
			doublereal Uoperating=VyMAX; // U_inf - �� �������� ����������� ������ �� �������������. (������ ������� �� ������������ �������� - VyMAX)
			while (f[0].neighbors_for_the_internal_node[T_SIDE][0][iBT]<f[0].maxelm) {
				volume3D(iBT, f[0].nvtx, f[0].pa, dx, dy, dz);
			    displacement_thickness[ilengthcounter]+=dz*(Uoperating-f[0].potent[VELOCITY_Y_COMPONENT][iBT]);
                iBT=f[0].neighbors_for_the_internal_node[T_SIDE][0][iBT]; // ��������� �� �������� ��������������� � ���������.
			}
			displacement_thickness[ilengthcounter]/=Uoperating; // ������� ����������.

			// ������������� ����������� ����.
			while (f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]<f[0].maxelm) iBT=f[0].neighbors_for_the_internal_node[B_SIDE][0][iBT]; // ��������� � ������ ��������.
			doublereal TempB, TempT;
			iBT=f[0].ptr[iBT];
			TempB=t.potent[t.neighbors_for_the_internal_node[B_SIDE][0][iBT]];
			TempT=t.potent[iBT];
			center_cord3D(iBT, t.nvtx, t.pa, p, B_SIDE); // ���������� ��������� ������ ��.
			volume3D(iBT, t.nvtx, t.pa, dx, dy, dz);
			zB=p.z-0.5*dz;
			zT=p.z;

			// ��������� ��������� ����� ���������:
			local_Nusselt_number[ilengthcounter]=(t.prop_b[LAM][t.neighbors_for_the_internal_node[B_SIDE][0][iBT]-t.maxelm]*((Twall-TempT)/(0.5*dz))*yposition[ilengthcounter])/((Twall-Tinf)*lam_avg);

			// ����������� ���������� �������������� ������������ ����.
			doublereal Temp_critical=Tinf+deltaTscal*(Twall-Tinf);
			while ((t.neighbors_for_the_internal_node[T_SIDE][0][iBT]<t.maxelm)&&(t.potent[iBT]>Temp_critical)) {
				iBT=t.neighbors_for_the_internal_node[T_SIDE][0][iBT]; // ��������� �� �������� ��������������� � ���������.
				TempB=t.potent[t.neighbors_for_the_internal_node[B_SIDE][0][iBT]];
			    TempT=t.potent[iBT];
				zB=zT;
                center_cord3D(iBT, t.nvtx, t.pa, p, T_SIDE); // ���������� ��������� ������ ��.
				zT=p.z;
			}
			a=(TempB-TempT)/(zB-zT);
			b=(TempT*zB-TempB*zT)/(zB-zT);
			deltaT[ilengthcounter]=(Temp_critical-b)/a; // ������� �������������� ������������ ����.

            ilengthcounter++;
		}
		iP=f[0].neighbors_for_the_internal_node[N_SIDE][0][iP];
	}

	// ������� ������������ ����������� ������ �� ������:
	// �������� ��������� �� ������� ��������.
	avg_Cx=0.0;
	doublereal slen=0.0;
	for (integer i=0; i<(iclength-1); i++) {
		slen+=(yposition[i+1]-yposition[i]);
		avg_Cx+=0.5*(yposition[i+1]-yposition[i])*(Cx[i+1]+Cx[i]);
	}
	avg_Cx=avg_Cx/slen; // ������� ������������ ����������� ������ �� ������.

	// ������ �������������� ���������� � ��������� ���� blasius_1908.txt
	FILE *fpblas=NULL; // ���� � ������� ����� ������������ ���������� � ������ ��������.
	
#ifdef MINGW_COMPILLER
	int err_blas = 0;
	fpblas=fopen64("blasius_1908.txt", "w");
	if (fpblas == NULL) err_blas = 1;
#else
	errno_t err_blas = 0;
	err_blas = fopen_s(&fpblas, "blasius_1908.txt", "w");
#endif

	if ((err_blas) != 0) {
		 printf("Create File blasius_1908.txt Error\n");
         //getchar();
		 system("pause");
         //exit(0);
     }
	 else {

		 if (fpblas != NULL) {
			 // ������ ���������� � ��������� ����:
			 fprintf(fpblas, "Laminar Flow and Heat Transfer over a Flat Plate.\n\n");
			 // �������� ����������:
			 fprintf(fpblas, "Average material properties: \n");
			 fprintf(fpblas, "Density= %+.16f kg/m!3\n", rho_avg);
			 fprintf(fpblas, "Dynamic_Viscosity= %+.16f Pa*s\n", mu_avg);
			 fprintf(fpblas, "Thermal_conductivity= %+.16f W/(m*K)\n", lam_avg);
			 fprintf(fpblas, "Specific_Heat= %+.16f J/(kg*K)\n\n", cp_avg);
			 // �������������� �������:
			 fprintf(fpblas, "Geometric_parametrs: \n");
			 fprintf(fpblas, "length_plate= %+.16f m\n\n", length_plate);
			 // ��������� �������:
			 fprintf(fpblas, "Boundary_conditions: \n");
			 fprintf(fpblas, "Inlet_fluid_velocity= %+.16f m/s\n", U_inf);
			 fprintf(fpblas, "Inlet_fluid_temperature= %+.16f oC\n", Tinf);
			 fprintf(fpblas, "Plate_Wall_temperature= %+.16f oC\n\n", Twall);
			 // ��������� ��������� ���������:
			 fprintf(fpblas, "Other_scalar_parametrs: \n");
			 fprintf(fpblas, "Critical_Renolds_number= %+.16f\n", 1e5);
			 fprintf(fpblas, "Average_Renolds_number= %+.16f\n", avg_Re_number);
			 fprintf(fpblas, "Average_Prandtl_number= %+.16f\n", avg_Pr_number);
			 fprintf(fpblas, "Average_Peclet_number= %+.16f\n", avg_Pe_number);
			 fprintf(fpblas, "Local_average_coefficient_of_friction_at_the_wall= %+.16f \n\n", avg_Cx);
			 // XY plot:
			 fprintf(fpblas, "XY Plots: \n");
			 fprintf(fpblas, "ypos - y_position m; \n");
			 fprintf(fpblas, "delta - 95%%_Boundary_layer_thickness m; \n");
			 fprintf(fpblas, "Cx - Wall_skin_friction_distribution; \n");
			 fprintf(fpblas, "delta* - displacement_thickness m; \n");
			 fprintf(fpblas, "deltaT - 5%%_Thermal_boundary_layer_distribution m; \n");
			 fprintf(fpblas, "local_Nusselt_number - Wall_Nusselt_number_distribution; \n\n");
			 // ���������� ���� �������� ������ ��� ���������� ��������:
			 fprintf(fpblas, "ypos	delta	Cx	delta*	deltaT	local_Nusselt_number\n");
			 for (integer i = 0; i < iclength; i++) {
				 fprintf(fpblas, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n", yposition[i], delta[i], Cx[i], displacement_thickness[i], deltaT[i], local_Nusselt_number[i]);
			 }

			 fclose(fpblas); // �������� ����� ��� ������ ��������������� � ������ �������� 1908 ����.
		 }
	 }


	// ������������ ����������� ������:

	// ����������������� ��������������
	delete[] delta; 
	delete[] yposition;
	delete[] Cx;
	delete[] displacement_thickness;
	// �������� ��������������
	delete[] deltaT;
	delete[] local_Nusselt_number;

} // boundarylayer_info

#endif
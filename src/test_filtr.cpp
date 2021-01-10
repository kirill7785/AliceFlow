// ���� test_filtr.cpp �������� �������� ���
// ��� ��������� ������ �������������� ���� ���������, � 
// ����� ����������� ������� �� ����: 
// 1. ��������� ������ �������������� ���� ��������,
// 2. ����� (������ ������ ��������) ��� ������ �������������� ���� ��������,
// 3. ���� ����� ������ � ������������� ������ (������������ ��� Selective Smagorinsky SS),
// 4. ���������� ������� ��������� ���������� ��� ������ �������������� ���� ���������.
// begin 19 ������ 2012 ����.

#ifndef MY_TEST_FILTR_CPP
#define MY_TEST_FILTR_CPP 1

// ������������ ������������ �� ����� ��������� ��������.
// 1) ������� �� ������, 2) Trapezoidal Filtr, 3) Simpson filtr.

// ����� ������� ������ ������ ������� �� ������ ��������.
// layer TOP
// 1	1	1
// 1	1	1
// 1	1	1
// layer center
// 1	1	1
// 1	1	1
// 1	1	1
// layer Bottom
// 1	1	1
// 1	1	1
// 1	1	1
// ��������� 1.0/27.0
#define VOLUME_AVERAGE_FILTR 0
// ������� Trapezoidal � SIMPSON �������� �� �������� (� ������) ������ ��� ����������� �����.
// Trapezoidal Filtr
// layer TOP
// 1	2	1
// 2	4	2
// 1	2	1
// layer center
// 2	4	2
// 4	8	4
// 2	4	2
// layer Bottom
// 1	2	1
// 2	4	2
// 1	2	1
// ��������� 1.0/64.0.
#define TRAPEZOIDAL_FILTR 1
// Simpson filtr
// layer TOP
// 1	4	1
// 4	16	4
// 1	4	1
// layer center
// 4	16	4
// 16	64	16
// 4	16	4
// layer Bottom
// 1	4	1
// 4	16	4
// 1	4	1
// ��������� 1.0/216.0.
#define SIMPSON_FILTR 2

// ������������ ����� flattener - ������������ �������
// ������������� � ����� ������ �� �������� ������ ����������� ��� 
// ������� �������� ��������� �� ����������� ��� �������� ����� ���������
// ����� ���� �����.
// ���������� ��� ������ ����� ������ ��������� ���������� ������� ��-�� ��� ����������������.
// �������� ! �� ��������� ����������� ����� ������ ����� �� �������� ������������ ������ potent.
// ������� �������������� ��� ������������ ���������� ����� ����������� ��������� �����.
// � ����� �� ����������� ��������� ���� ���. ������ ��� ������ ������ ��� ������������ ���� ��� �� �������� �� ��������.
// ���������������� ���� ����� ��������� ��������� ����� ����� ���������� ����������� (��� ��������� ��������).
void flattener(doublereal* &potent, integer maxelm, integer maxbound, int*** neighbors_for_the_internal_node,
	           int** nvtx, TOCHKA* pa, BOUND* border_neighbor) {

	doublereal* potentcopy = nullptr;
	potentcopy=new doublereal[maxelm + maxbound]; // ������� ������ ����������� ������������ ����.
	for (integer i=0; i<maxelm+maxbound; i++) {
		potentcopy[i]=potent[i];
	}

	bool** avgonsosed=new bool*[maxelm];
	for (integer i=0; i<maxelm; i++) avgonsosed[i]=new bool[6];
	for (integer i=0; i<maxelm; i++) for (integer j=0; j<6; j++) avgonsosed[i][j]=false;

	for (integer i=0; i<maxelm+maxbound; i++) {
		potent[i]=0.0; // ���������.
	}

	for (integer iP=0; iP<maxelm; iP++) {

		  doublereal dx, dy, dz;

		  
		  dx = 0.0; dy = 0.0; dz = 0.0; // ���������� ����� ��������
		
		  volume3D(iP, nvtx, pa, dx, dy, dz);

		// ���� �� ���� ���������� ����������� �������.
		  bool bE, bW, bN, bS, bT, bB;
		  bE=false; bW=false; bN=false; bS=false; bT=false; bB=false;

		  integer iE, iN, iT, iW, iS, iB; // ������ �������� ����������� �������
	      iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

          if (iE>=maxelm) bE=true; // ���� true �� ���� �������� ���������.
	      if (iN>=maxelm) bN=true;
	      if (iT>=maxelm) bT=true;
          if (iW>=maxelm) bW=true;
	      if (iS>=maxelm) bS=true;
	      if (iB>=maxelm) bB=true;

		  if (!bE) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		      volume3D(iE, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dx-dx1)>1e-20) avgonsosed[iP][E_SIDE]=true;
		  }

		  if (!bW) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		      volume3D(iW, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dx-dx1)>1e-20) avgonsosed[iP][W_SIDE]=true;
		  }

		  if (!bN) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		      volume3D(iN, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dy-dy1)>1e-20) avgonsosed[iP][N_SIDE]=true;
		  }

		  if (!bS) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		      volume3D(iS, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dy-dy1)>1e-20) avgonsosed[iP][S_SIDE]=true;
		  }

		  if (!bT) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		      volume3D(iT, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dz-dz1)>1e-20) avgonsosed[iP][T_SIDE]=true;
		  }

		  if (!bB) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		      volume3D(iB, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dz-dz1)>1e-20) avgonsosed[iP][B_SIDE]=true;
		  }
	}

	for (integer iP=0; iP<maxelm; iP++) {

		// ���� �� ���� ���������� ����������� �������.
		  bool bE, bW, bN, bS, bT, bB;
		  bE=false; bW=false; bN=false; bS=false; bT=false; bB=false;

		  integer iE, iN, iT, iW, iS, iB; // ������ �������� ����������� �������
	      iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

          if (iE>=maxelm) bE=true; // ���� true �� ���� �������� ���������.
	      if (iN>=maxelm) bN=true;
	      if (iT>=maxelm) bT=true;
          if (iW>=maxelm) bW=true;
	      if (iS>=maxelm) bS=true;
	      if (iB>=maxelm) bB=true;

		  doublereal dx, dy, dz;

		  dx=0.0; dy=0.0; dz=0.0; // ���������� ����� ��������
		  volume3D(iP, nvtx, pa, dx, dy, dz);

		if (!(avgonsosed[iP][E_SIDE]||avgonsosed[iP][W_SIDE]||avgonsosed[iP][N_SIDE]||avgonsosed[iP][S_SIDE]||avgonsosed[iP][T_SIDE]||avgonsosed[iP][B_SIDE])) {
			potent[iP]=potentcopy[iP]; // ��������� �������� ��� ��������� (��� ����������� ��� ����� ����������� �����).
		}
		else {
			// ������� �������������� - ������������� ������ ��������������� ��������.

			doublereal rn=0.0;
			if (avgonsosed[iP][E_SIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		        volume3D(iE, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iE]/(dx1*dy1*dz1));

				rn+=1.0;
			}

            if (avgonsosed[iP][W_SIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		        volume3D(iW, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iW]/(dx1*dy1*dz1));

				rn+=1.0;
			}

			if (avgonsosed[iP][N_SIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		        volume3D(iN, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iN]/(dx1*dy1*dz1));

				rn+=1.0;
			}

            if (avgonsosed[iP][S_SIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		        volume3D(iS, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iS]/(dx1*dy1*dz1));

				rn+=1.0;
			}

			if (avgonsosed[iP][T_SIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		        volume3D(iT, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iT]/(dx1*dy1*dz1));

				rn+=1.0;
			}

            if (avgonsosed[iP][B_SIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // ���������� ����� ��������
		        volume3D(iB, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iB]/(dx1*dy1*dz1));

				rn+=1.0;
			}

			potent[iP]/=rn; // �������� ������� ��������������.
		}
	}

	// ��� ��������� �������� �� ������� ����� ��������� ����� ������������ ������������.
	// ����� ��������� ������������ ������������ ��� 
	// ���������� ������������� �������� � ��������� ����������� �������.
	for (integer iG=0; iG<maxbound; iG++) {
		// ���� �� ���� ��������� ����������� �������.

		doublereal dx=0.0, dy=0.0, dz=0.0;// ����� �������� ������������ ������
		// ���������� �������� �������� ������������ ������:
	    volume3D(border_neighbor[iG].iI, nvtx, pa, dx, dy, dz);
		TOCHKA pp,pb,pbb;

		// ��������� ���������� �������
		switch (border_neighbor[iG].Norm) {
		  case E_SIDE: // ������� ������� W
			       // ������������ ������������.
			      
				   center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,W_SIDE);
				   center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,W_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb,E_SIDE);

				   potent[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent[neighbors_for_the_internal_node[E_SIDE][0][border_neighbor[iG].iII]], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x-0.5*dx);
				   
			       break;
		  case W_SIDE: // ������� ������� E
			       // ������������ ������������.

				   
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,E_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,EE_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb,W_SIDE);
					
			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent[neighbors_for_the_internal_node[W_SIDE][0][border_neighbor[iG].iII]], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x+0.5*dx);

			       break;
		  case N_SIDE: // ������� ������� S
			       // ������������ ������������.
			       
			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,S_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,SS_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb,N_SIDE);

			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent[neighbors_for_the_internal_node[N_SIDE][0][border_neighbor[iG].iII]], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y-0.5*dy);
			  
			       break;
		  case S_SIDE:// ������� ������� N
			       // ������������ ������������.

			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,N_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,NN_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb,S_SIDE);

			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent[neighbors_for_the_internal_node[S_SIDE][0][border_neighbor[iG].iII]], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y+0.5*dy);
			  
			       break;
		  case T_SIDE: // ������� ������� B
			       // ������������ ������������.

		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,B_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,BB_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb,T_SIDE);


			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent[neighbors_for_the_internal_node[T_SIDE][0][border_neighbor[iG].iII]], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z-0.5*dz);

			       break;
		  case B_SIDE: // ������� ������� T
			       // ������������ ������������.
			  
			      
		          center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,T_SIDE);
		          center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,TT_SIDE);
			      center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb,B_SIDE);
					
			      potent[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent[neighbors_for_the_internal_node[B_SIDE][0][border_neighbor[iG].iII]], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z+0.5*dz); 
				  
			       break;
		}

	}

	if (potentcopy != nullptr) {
		delete[] potentcopy;
	}

	for (integer i = 0; i < maxelm; i++) {
		delete[] avgonsosed[i];
	}
	delete[] avgonsosed;

} // flattener

// �� ������ 27 �������� ������ ������� ��, � ������� ������������ ������������
// ��������� �������� ������� � ����� � ������������ pzvezda.
/*
doublereal rinterpolFinposition(TOCHKA* c27, doublereal* potent27, TOCHKA pzvezda) {
	const integer dirP=0;
	const integer dirE=1;

} // rinterpolFinposition
*/

// ������ ��������� ������� ������������ � ����� ���������� ������ ����
// � ������� double_average_potent
void my_additional_to_calc_average(bool bflag, integer iNODE, doublereal rmultiplyer,
	                               int** nvtx, TOCHKA* pa, doublereal* potent_in,
								   doublereal &rsvol, doublereal &rsum, doublereal &rcol, 
								   doublereal &rcolnumber) {

	// rmultiplyer - ������� ����������� ��������� � ����� �������.

	if (!bflag&&(iNODE!=-1)) {
		doublereal dx, dy, dz;
        dx=0.0; dy=0.0; dz=0.0; // ���������� ����� ��������
		volume3D(iNODE, nvtx, pa, dx, dy, dz);
		//rsvol+=dx*dy*dz;
		rsvol+=rmultiplyer*dx*dy*dz;
		//->//rsvol=1.0;
		rsum+=rmultiplyer*potent_in[iNODE]*dx*dy*dz;
		//-->//rsum+=rmultiplyer*potent_in[iNODE];
		rcol+=rmultiplyer;
		rcolnumber+=1.0;
    }

} // my_additional_to_calc_average

// ��������� ������ �������������� ���� potent_test_filtr=test_filtr(potent_in),
// ����� ���������� ��������� �������. ��������� potent_in ����� �������� ����� ��������:
// ���������� ��������, ������������ � ��.
// 16 ��� 2012 ���� ������� ������������ ������������ ������� ������� �� �������.
// �������� ������ �� ����������������� ������������� ������������� �����.
// �� �������� �� ���� �����.
void double_average_potent(doublereal* potent_in, doublereal* &potent_test_filtr, 
	integer maxelm, integer maxbound, int*** neighbors_for_the_internal_node,
							 int** nvtx, TOCHKA* pa, doublereal* &delta_test_filtr,
							 integer itype_filtr, BOUND* border_neighbor, integer iquadraticinterpolboud) {

	// iquadraticinterpolboud - ������� ������������ ��������� ��� ����������� ������������� �������� ������� ������� �� �������.
	// ��������� �������� 0 - ������ �������� �������� �� ���������� ����������� ��  �� �������. 2 - ������������ ������������.

	// potent_in - �� ���� ������� ���������� �������� ������� �������� � ���������� 
	// ������� ��������� �����-������. ��� ��� �� ���� ��� ������������� ��������
	// Box �������� ��� ������� �������� ��� ����������� �����.
	// potent_test_filtr - �������� ���������� �������� ���������� �� potent_in ����������� ���������
	// Box ������� �������� � 27 ����������� ������� ���������� ������ ����������� ����� ������ � ������.

    // delta_test_filtr ������ ��������� ������� ��� ������ ���������� ��
	// ������ ������ (���������� ������������ 27 ����������� �������). 
	// �������� �������� �� 0 .. maxelm-1. ������ ���� � ������� ����� ������������ �������� nullptr.

	// ���������� itype_filtr - �������� �� ��� �������: 0 - ������� ������� �� ������, 
	// 1 - ������ �� ������ ������� ��������, 2 - ������ �� ������ ������� ��������.

	// �������������.
	for (integer i=0; i<maxelm+maxbound; i++) potent_test_filtr[i]=0.0;

	// ���� �� ���� ���������� ����������� �������.
	// � ����� ������������ �������� ���������� (�������� ������).
	for (integer iP=0; iP<maxelm; iP++) {

		 // ��������. ��� ������� ����������� �� ����� ������������ �� ����� 27 (������� ��� ������) �������� (�� ���������) �������.
		 // ������ ��� ��������� ������� �������� ������� �������������� ������ ��� ������� �������������� ���������� �� ����� �� �����
		 // ���� �������� ����������� ������� ���������� ������ ����� ���������� ����� �� 26 ��������� ������� ������������� ������������.
		 // ����������� ����������� ������� � �������� "������������" ��������� �������� �� ���� ��������. ����� ��� ������ �������� ��������� 
		 // ����� ������� �������� ������� ����� �������: ���������� ����� � ������� � �� ����� �����. �.�. ��� ������� �� ������ ��������.
		 // ��� ����� ������������ - ��������� ��������. ��� ����� ������������ ������ �������� �������� - ��� �������� ���������� � ������ �� 
		 // ��������������� ��� � �������� ���������������� �� ���� ����������� �����.

          bool bE, bW, bN, bS, bT, bB;
		  bE=false; bW=false; bN=false; bS=false; bT=false; bB=false;

		  integer iE, iN, iT, iW, iS, iB; // ������ �������� ����������� �������
	      iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

          if (iE>=maxelm) bE=true; // ���� true �� ���� �������� ���������.
	      if (iN>=maxelm) bN=true;
	      if (iT>=maxelm) bT=true;
          if (iW>=maxelm) bW=true;
	      if (iS>=maxelm) bS=true;
	      if (iB>=maxelm) bB=true;

		  integer iEN=-1, iES=-1, iWS=-1, iWN=-1, iNT=-1, iNB=-1, iET=-1, iEB=-1, iST=-1, iSB=-1, iWT=-1, iWB=-1; // ������ �������� ����������� ������� (12 ��).
          bool bEN, bES, bWS, bWN, bNT, bNB, bET, bEB, bST, bSB, bWT, bWB; 
		  // false �������� ��� �� ����������, � true ��� �� ��������� ���� ��� ������ �� ����������.
		  bEN=false; bES=false; bWS=false; bWN=false; bNT=false; bNB=false;
		  bET=false; bEB=false; bST=false; bSB=false; bWT=false; bWB=false;

		  if (!bE && !bN) {
			  // ���������� ��.
			  iEN=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]];
			  if (iEN>=maxelm) bEN=true; // ��������� ��.
		  } else {
			  bEN=true;
			  iEN=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bE && !bS) {
			  // ���������� ��.
			  iES=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]];
			  if (iES>=maxelm) bES=true; // ��������� ��.
		  } else {
			  bES=true;
			  iES=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bW && !bS) {
			  // ���������� ��.
			  iWS=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][iP]];
			  if (iWS>=maxelm) bWS=true; // ��������� ��.
		  } else {
			  bWS=true;
			  iWS=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bW && !bN) {
			  // ���������� ��.
			  iWN=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][iP]];
			  if (iWN>=maxelm) bWN=true; // ��������� ��.
		  } else {
			  bWN=true;
			  iWN=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bN && !bT) {
			  // ���������� ��.
			  iNT=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]];
			  if (iNT>=maxelm) bNT=true; // ��������� ��.
		  } else {
			  bNT=true;
			  iNT=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bN && !bB) {
			  // ���������� ��.
			  iNB=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]];
			  if (iNB>=maxelm) bNB=true; // ��������� ��.
		  } else {
			  bNB=true;
			  iNB=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bE && !bT) {
			  // ���������� ��.
			  iET=neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]];
			  if (iET>=maxelm) bET=true; // ��������� ��.
		  } else {
			  bET=true;
			  iET=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bE && !bB) {
			  // ���������� ��.
			  iEB=neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]];
			  if (iEB>=maxelm) bEB=true; // ��������� ��.
		  } else {
			  bEB=true;
			  iEB=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bS && !bT) {
			  // ���������� ��.
			  iST=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]];
			  if (iST>=maxelm) bST=true; // ��������� ��.
		  } else {
			  bST=true;
			  iST=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bS && !bB) {
			  // ���������� ��.
			  iSB=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]];
			  if (iSB>=maxelm) bSB=true; // ��������� ��.
		  } else {
			  bSB=true;
			  iSB=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bW && !bT) {
			  // ���������� ��.
			  iWT=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]];
			  if (iWT>=maxelm) bWT=true; // ��������� ��.
		  } else {
			  bWT=true;
			  iWT=maxelm+maxbound; // KO �� ����������.
		  }

		  if (!bW && !bB) {
			  // ���������� ��.
			  iWB=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]];
			  if (iWB>=maxelm) bWB=true; // ��������� ��.
		  } else {
			  bWB=true;
			  iWB=maxelm+maxbound; // KO �� ����������.
		  }

		  integer iTNE=-1, iTSE=-1, iTNW=-1, iTSW=-1, iBNE=-1, iBSE=-1, iBNW=-1, iBSW=-1; // ������� ����� ���� (8 ��).
		  bool bTNE, bTSE, bTNW, bTSW, bBNE, bBSE, bBNW, bBSW;
          bTNE=true; bTSE=true; bTNW=true; bTSW=true; bBNE=true; bBSE=true; bBNW=true; bBSW=true; // �� ��������� ������ �� �� ����������,
		  // �� ���� ��������� ���� ��� �� ����������.

		  if (!bEN) {
			  if (bTNE) {
				  // ������������������ ������� ����� �����.
				  iTNE=neighbors_for_the_internal_node[T_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]]];
			      bTNE=false;
			      if (iTNE>=maxelm) bTNE=true;
			  }
			  if (bBNE) {
				  // ������������������ ������� ����� �����.
				  iBNE=neighbors_for_the_internal_node[B_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]]];
				  bBNE=false;
			      if (iBNE>=maxelm) bBNE=true;
			  }
		  }

          if (!bES) {
			  if (bTSE) {
				  // ������������������ ������� ����� �����.
				  iTSE=neighbors_for_the_internal_node[T_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]]];
			      bTSE=false;
			      if (iTSE>=maxelm) bTSE=true;
			  }
			  if (bBSE) {
				  // ������������������ ������� ����� �����.
				  iBSE=neighbors_for_the_internal_node[B_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]]];
				  bBSE=false;
			      if (iBSE>=maxelm) bBSE=true;
			  }
		  }

		  if (!bWN) {
			  if (bTNW) {
				  // ������������������ ������� ����� �����.
				  iTNW=neighbors_for_the_internal_node[T_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][iP]]];
			      bTNW=false;
			      if (iTNW>=maxelm) bTNW=true;
			  }
			  if (bBNW) {
				  // ������������������ ������� ����� �����.
				  iBNW=neighbors_for_the_internal_node[B_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][iP]]];
				  bBNW=false;
			      if (iBNW>=maxelm) bBNW=true;
			  }
		  }

          if (!bWS) {
			  if (bTSW) {
				  // ������������������ ������� ����� �����.
				  iTSW=neighbors_for_the_internal_node[T_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][iP]]];
			      bTSW=false;
			      if (iTSW>=maxelm) bTSW=true;
			  }
			  if (bBSW) {
				  // ������������������ ������� ����� �����.
				  iBSW=neighbors_for_the_internal_node[B_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][iP]]];
				  bBSW=false;
			      if (iBSW>=maxelm) bBSW=true;
			  }
		  }

		  if (!bNT) {
			  if (bTNE) {
				  // ������������������ ������� ����� �����.
				  iTNE=neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]]];
			      bTNE=false;
			      if (iTNE>=maxelm) bTNE=true;
			  }
			  if (bTNW) {
				  // ������������������ ������� ����� �����.
				  iTNW=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]]];
				  bTNW=false;
			      if (iTNW>=maxelm) bTNW=true;
			  }
		  }

		  if (!bNB) {
			  if (bBNE) {
				  // ������������������ ������� ����� �����.
				  iBNE=neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
			      bBNE=false;
			      if (iBNE>=maxelm) bBNE=true;
			  }
			  if (bBNW) {
				  // ������������������ ������� ����� �����.
				  iBNW=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
				  bBNW=false;
			      if (iBNW>=maxelm) bBNW=true;
			  }
		  }

		  if (!bET) {
			  if (bTNE) {
				  // ������������������ ������� ����� �����.
				  iTNE=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]]];
			      bTNE=false;
			      if (iTNE>=maxelm) bTNE=true;
			  }
			  if (bTSE) {
				  // ������������������ ������� ����� �����.
				  iTSE=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][iP]]];
				  bTSE=false;
			      if (iTSE>=maxelm) bTSE=true;
			  }
		  }

		  if (!bEB) {
			  if (bBNE) {
				  // ������������������ ������� ����� �����.
				  iBNE=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
			      bBNE=false;
			      if (iBNE>=maxelm) bBNE=true;
			  }
			  if (bBSE) {
				  // ������������������ ������� ����� �����.
				  iBSE=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
				  bBSE=false;
			      if (iBSE>=maxelm) bBSE=true;
			  }
		  }

          if (!bST) {
			  if (bTSE) {
				  // ������������������ ������� ����� �����.
				  iTSE=neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]]];
			      bTSE=false;
			      if (iTSE>=maxelm) bTSE=true;
			  }
			  if (bTSW) {
				  // ������������������ ������� ����� �����.
				  iTSW=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]]];
				  bTSW=false;
			      if (iTSW>=maxelm) bTSW=true;
			  }
		  }

		  if (!bSB) {
			  if (bBSE) {
				  // ������������������ ������� ����� �����.
				  iBSE=neighbors_for_the_internal_node[E_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
			      bBSE=false;
			      if (iBSE>=maxelm) bBSE=true;
			  }
			  if (bBSW) {
				  // ������������������ ������� ����� �����.
				  iBSW=neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
				  bBSW=false;
			      if (iBSW>=maxelm) bBSW=true;
			  }
		  }

		  if (!bWT) {
			  if (bTNW) {
				  // ������������������ ������� ����� �����.
				  iTNW=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]]];
			      bTNW=false;
			      if (iTNW>=maxelm) bTNW=true;
			  }
			  if (bTSW) {
				  // ������������������ ������� ����� �����.
				  iTSW=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[T_SIDE][0][iP]]];
				  bTSW=false;
			      if (iTSW>=maxelm) bTSW=true;
			  }
		  }

		  if (!bWB) {
			  if (bBNW) {
				  // ������������������ ������� ����� �����.
				  iBNW=neighbors_for_the_internal_node[N_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
			      bBNW=false;
			      if (iBNW>=maxelm) bBNW=true;
			  }
			  if (bBSW) {
				  // ������������������ ������� ����� �����.
				  iBSW=neighbors_for_the_internal_node[S_SIDE][0][neighbors_for_the_internal_node[W_SIDE][0][neighbors_for_the_internal_node[B_SIDE][0][iP]]];
				  bBSW=false;
			      if (iBSW>=maxelm) bBSW=true;
			  }
		  }

          // �� ��� ��������� 26 ���������� �� ������ ��������, ����� ���������� � ����������.
		  // potent_in

		  doublereal rsum=0.0, rsvol=0.0, rcol=0.0, rcolnumber=0.0;
		  doublereal dx, dy, dz;

		  dx=0.0; dy=0.0; dz=0.0; // ���������� ����� ��������
		  volume3D(iP, nvtx, pa, dx, dy, dz);
		  		  
		  //rsvol+=dx*dy*dz;	
		  switch (itype_filtr) {
		     case VOLUME_AVERAGE_FILTR: rcol+=1.0; rsum+=1.0*potent_in[iP]*dx*dy*dz; rsvol+=dx*dy*dz; break;
             case TRAPEZOIDAL_FILTR: rcol+=8.0; rsum+=8.0*potent_in[iP]*dx*dy*dz; rsvol+=8.0*dx*dy*dz; break;
	         case SIMPSON_FILTR: rcol+=64.0; rsum+=64.0*potent_in[iP]*dx*dy*dz; rsvol+=64.0*dx*dy*dz; break;
		  }
		  
		  /*
		  rsvol=1.0;
		  switch (itype_filtr) {
		     case VOLUME_AVERAGE_FILTR: rcol+=1.0; rsum+=1.0*potent_in[iP]; break;
             case TRAPEZOIDAL_FILTR: rcol+=8.0; rsum+=8.0*potent_in[iP]; break;
	         case SIMPSON_FILTR: rcol+=64.0; rsum+=64.0*potent_in[iP]; break;
		  }
		  */
		  rcolnumber+=1.0;

		  /*
		  // ���� � ����������� ��� 
		  if (!bE) {
              dx=0.0; dy=0.0; dz=0.0; // ���������� ����� ��������
		      volume3D(iE, nvtx, pa, dx, dy, dz);
		      rsvol+=dx*dy*dz;
		      rsum+=potent_in[iE]*dx*dy*dz; 
		  }
		  ������ �� ����� �������
		   my_additional_to_calc_average(bE, iE, nvtx, pa, potent_in, rsvol, rsum);
		   // ��� �������� ������� ����������.
		  */

		  switch (itype_filtr) {
			  case VOLUME_AVERAGE_FILTR:

				  // ������� ������� �� ������.
				  
		           // ������� ��������� ������. (6��)
		           my_additional_to_calc_average(bE, iE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bW, iW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bN, iN, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bS, iS, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bT, iT, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bB, iB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // ������ ������ ������� ������ (12��).
		           my_additional_to_calc_average(bEN, iEN, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bES, iES, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWS, iWS, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWN, iWN, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNT, iNT, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNB, iNB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bET, iET, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bEB, iEB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bST, iST, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bSB, iSB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWT, iWT, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWB, iWB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // ������ ����� ������� ������ (8��).
		           my_additional_to_calc_average(bTNE, iTNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSE, iTSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTNW, iTNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSW, iTSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNE, iBNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSE, iBSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNW, iBNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSW, iBSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);  

				   break;

		      case TRAPEZOIDAL_FILTR: 
				  
				   // ������������� ������ ���������� �� ������� ��������.
				   // ������� ������������ �������� ��� ����������� �����.

				   // ������� ��������� ������. (6��)
		           my_additional_to_calc_average(bE, iE, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bW, iW, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bN, iN, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bS, iS, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bT, iT, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bB, iB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // ������ ������ ������� ������ (12��).
		           my_additional_to_calc_average(bEN, iEN, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bES, iES, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWS, iWS, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWN, iWN, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNT, iNT, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNB, iNB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bET, iET, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bEB, iEB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bST, iST, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bSB, iSB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWT, iWT, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWB, iWB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // ������ ����� ������� ������ (8��).
		           my_additional_to_calc_average(bTNE, iTNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSE, iTSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTNW, iTNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSW, iTSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNE, iBNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSE, iBSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNW, iBNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSW, iBSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);

				   break;

              case SIMPSON_FILTR:

				   // ������ ���������� �� ������� ��������.
				   // ������� ������������ �������� ��� ����������� �����.

				   // ������� ��������� ������. (6��)
		           my_additional_to_calc_average(bE, iE, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bW, iW, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bN, iN, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bS, iS, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bT, iT, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bB, iB, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // ������ ������ ������� ������ (12��).
		           my_additional_to_calc_average(bEN, iEN, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bES, iES, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWS, iWS, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWN, iWN, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNT, iNT, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNB, iNB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bET, iET, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bEB, iEB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bST, iST, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bSB, iSB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWT, iWT, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWB, iWB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // ������ ����� ������� ������ (8��).
		           my_additional_to_calc_average(bTNE, iTNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSE, iTSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTNW, iTNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSW, iTSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNE, iBNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSE, iBSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNW, iBNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSW, iBSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber); 

				   break;
		  }

		  if (delta_test_filtr!=nullptr) {
             delta_test_filtr[iP]=exp((1.0/3.0)*log(rsvol)); // ������ ��������� ������� � ������ iP.
		  }

		  // ���� ���� ������� ������� ���� ��� ������ �������.
		  //potent_test_filtr[iP]=(rsum/rcol)*(rcolnumber/rsvol); // ������� �� ������ ��������.
          potent_test_filtr[iP]=(rsum/rsvol); // ������� �� ������ ��������.
		  // �� ��� ������ �������� �������� ����� ��������� ������������ ������������ �� ����������� ����� � �������������.
		  //potent_test_filtr[iP]=(rsum/rcol);//*(rcolnumber); // ������ ������� ��� ����� ������.
	}

	
	if (iquadraticinterpolboud==2) {
		// ����� ��������� ������������ ������������ ��� 
	// ���������� ������������� �������� � ��������� ����������� �������.
	for (integer iG=0; iG<maxbound; iG++) {
		// ���� �� ���� ��������� ����������� �������.

		doublereal dx=0.0, dy=0.0, dz=0.0;// ����� �������� ������������ ������
		// ���������� �������� �������� ������������ ������:
	    volume3D(border_neighbor[iG].iI, nvtx, pa, dx, dy, dz);
		TOCHKA pp,pb,pbb;

		// ��������� ���������� �������
		switch (border_neighbor[iG].Norm) {
		  case E_SIDE: // ������� ������� W
			       // ������������ ������������.
			      
				   center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, W_SIDE);
				   center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, WW_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[E_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb, E_SIDE);

				   potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent_test_filtr[neighbors_for_the_internal_node[E_SIDE][0][border_neighbor[iG].iII]], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x-0.5*dx);
				   
			       break;
		  case W_SIDE: // ������� ������� E
			       // ������������ ������������.

				   
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, E_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, EE_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[W_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb, W_SIDE);
					
			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent_test_filtr[neighbors_for_the_internal_node[W_SIDE][0][border_neighbor[iG].iII]], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x+0.5*dx);

			       break;
		  case N_SIDE: // ������� ������� S
			       // ������������ ������������.
			       
			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, S_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, SS_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[N_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb, N_SIDE);

			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent_test_filtr[neighbors_for_the_internal_node[N_SIDE][0][border_neighbor[iG].iII]], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y-0.5*dy);
			  
			       break;
		  case S_SIDE:// ������� ������� N
			       // ������������ ������������.

			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, N_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,NN_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[S_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb, S_SIDE);

			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent_test_filtr[neighbors_for_the_internal_node[S_SIDE][0][border_neighbor[iG].iII]], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y+0.5*dy);
			  
			       break;
		  case T_SIDE: // ������� ������� B
			       // ������������ ������������.

		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, B_SIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, BB_SIDE);
			       center_cord3D(neighbors_for_the_internal_node[T_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb, T_SIDE);


			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent_test_filtr[neighbors_for_the_internal_node[T_SIDE][0][border_neighbor[iG].iII]], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z-0.5*dz);

			       break;
		  case B_SIDE: // ������� ������� T
			       // ������������ ������������.
			  
			      
		          center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, T_SIDE);
		          center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, TT_SIDE);
			      center_cord3D(neighbors_for_the_internal_node[B_SIDE][0][border_neighbor[iG].iII], nvtx, pa, pbb, B_SIDE);
					
			      potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent_test_filtr[neighbors_for_the_internal_node[B_SIDE][0][border_neighbor[iG].iII]], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z+0.5*dz); 
				  
			       break;
		}

	}
	}
	else {
		for (integer iG=0; iG<maxbound; iG++) {
		   // ���� �� ���� ��������� ����������� �������.
           potent_test_filtr[border_neighbor[iG].iB]=potent_test_filtr[border_neighbor[iG].iI];

	    }
	}


} // double_average_velocity

// ��� ������� ����������� ������������ ������ ���������� �������� true ���� �������� ���� ����� ������ � 
// ������������� ������ ������ ��� beta0==15 ��������. ������������ � ������ �������������� Selective Smagorinsky (SS).
void calc_selective_smagorinsky(FLOW &f, bool* &bfibeta, integer itype_filtr, doublereal beta0) {

	// �������� ������� ���� �����.
	// ������ �������������� Selective Smagorinsky ���������� �� ��������� ����� ������������ ������������ ��������
	// � ��� ������ ���� ���� ����� ������ � ���������� ������ ������ ���������� �������� beta0. �������� beta0 �����������
	// � ������������ � �������� ���������� ����� 15 ��������. ������ ���� �������������� ����� ������ ������������ ������������ 
	// �������� ��������� ������������ ��������� ��������������� ������ �������������.
	// ������� ���������: 
	// f - ��� ���������� � ����������������� ����������,
	// beta0 - ��������� �������� ���� ����� ������ � ���������� ������ � ������������� ������ �������������,
	// itype_filtr - ���� �� ��� �������� �� ����� , ��������, SIMPSON_FILTR - ������ ���� ��������.
	// bfibeta - ������������ ��� ���� ������ ���������� ����������� ������� ������ ��������, ��������� � ��� 
	// ����� �� ��������� ����� ������������ ������������ �������� ��� ���.

	// ��������������� ��������:
	// ������� ���������� ������������:
	// [vectora x vectorb]=(ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx);

	// bfibeta - ��������� �������� true (������������ ������������ �������� ����� ���������) � �������� false
	// ���� ������������ ������������ �������� ��������� �������. 
	// bfibeta ����� ����������� 0..maxelm-1 �.�. ����� �������� ��� ���� ���������� ����������� �������.

	// �������� �����:
	// 1. ����� ��� ���������� ��������.
	// 2. ���� ������� ����������� ��������� �������� �� �.1 ��������� ��� ���������� ����� �� ������� ���������� ������������.
	// 3. ������������� � ��������� ��������������� ���������� �����.
	// 4. ��������� ������ ���������� ������������ ����� �� ������������� ����� ��� ������ �� ����� ��������� ��������� ���������������� �������.
	// 5. ����� ������ ����� � ������ �������������� �����. 
	// 6. ��������� ������� ���� � �������� �������� ������� ��������.

	// VX - 0, VY - 1, VZ - 2

	
	// ���������� �����.
	const integer iCURLX=0;
	const integer iCURLY=1;
	const integer iCURLZ=2;

	doublereal** curl_component = nullptr;
	curl_component = new doublereal*[3];
	for (integer i=0; i<3; i++) curl_component[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		// ������������� .
		curl_component[iCURLX][i]=0.0;
		curl_component[iCURLY][i]=0.0;
		curl_component[iCURLZ][i]=0.0;
	}

	// ���������� �����.
    for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		curl_component[iCURLX][i]=f.potent[GRADYVZ][i]-f.potent[GRADZVY][i];
		curl_component[iCURLY][i]=f.potent[GRADZVX][i]-f.potent[GRADXVZ][i];
	    curl_component[iCURLZ][i]=f.potent[GRADXVY][i]-f.potent[GRADYVX][i];
	}

	// ��������� ����������� ������ ��� ��������������� ���������� �����.
	doublereal** curl_component_test_filtr = nullptr;
	curl_component_test_filtr = new doublereal*[3];
	for (integer i=0; i<3; i++) curl_component_test_filtr[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		// ������������� .
		curl_component_test_filtr[iCURLX][i]=0.0;
		curl_component_test_filtr[iCURLY][i]=0.0;
		curl_component_test_filtr[iCURLZ][i]=0.0;
	}

	doublereal* delta_test_filtr;
	delta_test_filtr=nullptr;

	// ���������� ������������� ��������� �����:
	// �������� ! ������ �� ����������� ��������� �������� �� ����� ����������� �� ���� ������������� ��������.
	// ��� ���������� �������� �� ������� ������� ����������� ������������ ������������ ������� ������� ������.
	double_average_potent(curl_component[iCURLX], curl_component_test_filtr[iCURLX], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
		                     f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(curl_component[iCURLY], curl_component_test_filtr[iCURLY], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(curl_component[iCURLZ], curl_component_test_filtr[iCURLZ], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	bfibeta=new bool[f.maxelm]; // ��������� ����������� ������ ��� ������������ ��������.
	for (integer i=0; i<f.maxelm; i++) bfibeta[i]=true; // �� ��������� ������������ �������� ����� ����� ���������.
	beta0=fabs(beta0);

	doublereal* curlxtest_curl = nullptr;
	curlxtest_curl = new doublereal[3];
	// ���������� �������� ������������ ������� �������:
	for (integer i=0; i<f.maxelm; i++) {
		// ��� ���� ���������� ����������� �������.
		doublereal module_Curl=0.0;
		// ������ �����.
		module_Curl=sqrt(curl_component[iCURLX][i]*curl_component[iCURLX][i]+curl_component[iCURLY][i]*curl_component[iCURLY][i]+curl_component[iCURLZ][i]*curl_component[iCURLZ][i]);

        doublereal module_Curl_test_filtr=0.0;
		// ������ �������������� �����.
		module_Curl_test_filtr=sqrt(curl_component_test_filtr[iCURLX][i]*curl_component_test_filtr[iCURLX][i]+curl_component_test_filtr[iCURLY][i]*curl_component_test_filtr[iCURLY][i]+curl_component_test_filtr[iCURLZ][i]*curl_component_test_filtr[iCURLZ][i]);

		
		curlxtest_curl[iCURLX]=curl_component[iCURLY][i]*curl_component_test_filtr[iCURLZ][i]-curl_component[iCURLZ][i]*curl_component_test_filtr[iCURLY][i];
        curlxtest_curl[iCURLY]=curl_component[iCURLZ][i]*curl_component_test_filtr[iCURLX][i]-curl_component[iCURLX][i]*curl_component_test_filtr[iCURLZ][i];
		curlxtest_curl[iCURLZ]=curl_component[iCURLX][i]*curl_component_test_filtr[iCURLY][i]-curl_component[iCURLY][i]*curl_component_test_filtr[iCURLX][i];
		doublereal module_curlxcurl=0.0;
        module_curlxcurl=sqrt(curlxtest_curl[iCURLX]*curlxtest_curl[iCURLX]+curlxtest_curl[iCURLY]*curlxtest_curl[iCURLY]+curlxtest_curl[iCURLZ]*curlxtest_curl[iCURLZ]);
		

		// ��� ���������� ���� ����� ��������� �������������� ��������,
		// ��������, ������� �� ���� � �.�. ������� ���� �������� ������������
		// � ����� ������ � ����� �� ����������. 
		try
		{
			if ((module_Curl<1e-30) || (module_Curl_test_filtr<1e-30)) {
				// �����-�� �������� ����� ������������ ���� ��������
				// ������������ �������� ������������ ������������ ��������.
                bfibeta[i]=false;
				// ��� �������� ����� �������� �� ������ �������� ����� ����� ����� ����� ����.
				//printf("warning little module Curl...\n");
				//getchar();
			}
			else {
				// �������� ��������� ������ ������ �� -1.0 �� +1.0.
			    doublereal beta=180.0*fabs(asin(module_curlxcurl/(module_Curl*module_Curl_test_filtr)))/M_PI;
			    if (beta>=beta0) {
				   // ���� ���������� � ������ ����� ��������� ������������ ��������.
				   bfibeta[i]=true;
			    }
			    else {
				   // ���� ������� ��������� � ������������ �������� 
				   // ��������� �������.
				   bfibeta[i]=false;
			    }
			}
		}
		catch (integer a)
		{
			printf("Selective Smagorinsky angle beta exeption...\n");
#if doubleintprecision == 1
			printf("Caught exception number:  %lld\n", a);
#else
			printf("Caught exception number:  %d\n", a);
#endif

			
			printf("Please, press any key to exit...\n");
			//getchar();
			system("pause");
			exit(0);
		}
		
	}

	if (curlxtest_curl != nullptr) {
		delete[] curlxtest_curl;
	}
	// ������������ ����������� ������.
	// ����������� �����.
	if (curl_component != nullptr) {
		for (integer i = 0; i < 3; i++) {
			if (curl_component[i] != nullptr) {
				delete[] curl_component[i];
			}
		}
		delete[] curl_component;
	}
	// ����������� ���������������� �����.
	if (curl_component_test_filtr != nullptr) {
		for (integer i = 0; i < 3; i++) {
			if (curl_component_test_filtr[i] != nullptr) {
				delete[] curl_component_test_filtr[i];
			}
		}
		delete[] curl_component_test_filtr;
	}

} // calc_selective_smagorinsky


// begin 28 ������ 2012 ����.
// ���������� ������ �������������� �������.
// Dynamic subgrid-scale model 1991 ���.
// ��� ������������ ����� ������ �������������� 1991 ���.
// ���� ������ ������� � ���, ��� ��� ���� ������������ �������������
// ������������� (������� �� �������� ����������� ��� � 1960 ����� ��� ����������������
// ������������ �� ������ 20x20x20) � ���������� �������� ��� ��������� Cs ������� � 
// ������� �� ���� ��������� Cg - ��������� �������. �������� ��������� ������� �����������
// �� ������ ���������� ������������ � ��������� ����� ������� �������� � �������� ��������.
// ��� ������� �������� ���������� ������� ������������ �������� ������ ������������ ������.
// ��. ����� [1] �.�. �� ������ � �������� ������������� ������������ ������� � ������������, 
// ���������, ����������� ��������� � ���������� �������. ������ 2009 ����.
// [2] ������ �.�., ��������� �.�. ������������� ������� ������ � �������� ������������ �������.
// ���������, 2008.
void my_Germano_model(FLOW &f, doublereal* &Cs2, integer itype_filtr) {
	// ������ ����� ���������� ������� ��������� �������������.
	// �������� �������� �������� ����� ���� � ������������� ��� ���� ����� ����������
	// ����� ��� �������� ������� �� ������� �� ������ ������ � ������� (�.�. � ���������������
	// ������� �� ������������� ����������� �������� �������).

	// �������� Cs2 - ���������� ������ ��� ���������� ����������� �������.
	// ����������� 0..f.maxelm-1.

	// ������ ���� ������� ������ ��� Cs2 ���������� �������������.

	// �������� ����� ��������� ������������ ��� ������� ��������� ��������
	// ������������� ��������.


	// VX 0 VY 1 VZ 2
	doublereal** speed_test_filtering=new doublereal*[3];
	for (integer i=0; i<3; i++) {
		speed_test_filtering[i]=new doublereal[f.maxelm+f.maxbound];
	}
	// �������������
	for (integer i=0; i<3; i++) {
		#pragma omp parallel for
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			speed_test_filtering[i][j]=0.0;  // �������������.
		}
	}

	doublereal* delta_test_filtr;
	delta_test_filtr=nullptr;


	// ���������� ������������� ��������� ��������.
	// �������� ! ������ �� ����������� ��������� �������� �� ����� ����������� �� ���� ������������� ��������.
	// ��� ���������� �������� �� ������� ������� ����������� ������������ ������������ ������� ������� ������.
	double_average_potent(f.potent[VELOCITY_X_COMPONENT], speed_test_filtering[VELOCITY_X_COMPONENT], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(f.potent[VELOCITY_Y_COMPONENT], speed_test_filtering[VELOCITY_Y_COMPONENT], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(f.potent[VELOCITY_Z_COMPONENT], speed_test_filtering[VELOCITY_Z_COMPONENT], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	

	// ��� ���������� ���������� �������� ����� ������� �� �������� ������������ ��������� ��������.
	const integer VXVX=0;
	const integer VXVY=1;
	const integer VXVZ=2;
	const integer VYVX=3;
	const integer VYVY=4;
	const integer VYVZ=5;
	const integer VZVX=6;
	const integer VZVY=7;
	const integer VZVZ=8;

	doublereal** speedxspeed_test_filtering = nullptr;
	speedxspeed_test_filtering = new doublereal*[9];
	for (integer i=0; i<9; i++) {
		speedxspeed_test_filtering[i]=new doublereal[f.maxelm+f.maxbound];
	}
	// �������������
	for (integer i=0; i<9; i++) {
		#pragma omp parallel for
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			speedxspeed_test_filtering[i][j]=0.0;  // �������������.
		}
	}

	doublereal* speedxspeed_buf = nullptr;
	speedxspeed_buf = new doublereal[f.maxelm + f.maxbound];
	// �������������
	#pragma omp parallel for
	for (integer j=0; j<f.maxelm+f.maxbound; j++) {
		speedxspeed_buf[j]=0.0;
	}

	// ���������� �������� ������������:
	integer** imarker=new integer*[9];
	for (integer i=0; i<9; i++) imarker[i]=new integer[2];
	imarker[VXVX][0]=VELOCITY_X_COMPONENT; imarker[VXVX][1]=VELOCITY_X_COMPONENT;
	imarker[VXVY][0]=VELOCITY_X_COMPONENT; imarker[VXVY][1]=VELOCITY_Y_COMPONENT;
	imarker[VXVZ][0]=VELOCITY_X_COMPONENT; imarker[VXVZ][1]=VELOCITY_Z_COMPONENT;
	imarker[VYVX][0]=VELOCITY_Y_COMPONENT; imarker[VYVX][1]=VELOCITY_X_COMPONENT;
	imarker[VYVY][0]=VELOCITY_Y_COMPONENT; imarker[VYVY][1]=VELOCITY_Y_COMPONENT;
	imarker[VYVZ][0]=VELOCITY_Y_COMPONENT; imarker[VYVZ][1]=VELOCITY_Z_COMPONENT;
	imarker[VZVX][0]=VELOCITY_Z_COMPONENT; imarker[VZVX][1]=VELOCITY_X_COMPONENT;
	imarker[VZVY][0]=VELOCITY_Z_COMPONENT; imarker[VZVY][1]=VELOCITY_Y_COMPONENT;
	imarker[VZVZ][0]=VELOCITY_Z_COMPONENT; imarker[VZVZ][1]=VELOCITY_Z_COMPONENT;

    for (integer i=0; i<9; i++) {
		integer fmr=imarker[i][0], smr=imarker[i][1]; // first marker, second marker

	    #pragma omp parallel for
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
		    speedxspeed_buf[j]=f.potent[fmr][j]*f.potent[smr][j];
	    }
		// ���������� ������������� �������� ������������ ��������� ��������.
	    // �������� ! ������ �� ����������� ��������� �������� �� ����� ����������� �� ���� ������������� ��������.
	    // ��� ���������� �������� �� ������� ������� ����������� ������������ ������������ ������� ������� ������.
	    double_average_potent(speedxspeed_buf, speedxspeed_test_filtering[i], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	}

	// ���������� ���������� ��������.
	doublereal** rLeonard_stress = nullptr;
	rLeonard_stress = new doublereal*[9];
	for (integer i=0; i<9; i++) {
		rLeonard_stress[i]=new doublereal[f.maxelm+f.maxbound];
	}
	// ���������������� ���������� ���������� �������� (1974):
	// ��. ������� (2.20) �� 43 �������� ����� �.�. ��� ��. [1].
	for (integer i=0; i<9; i++) {
		integer fmr=imarker[i][0], smr=imarker[i][1]; // first marker, second marker

		#pragma omp parallel for 
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			// ��� ����������� ����� �� ����� ���. 
			// ��� ����������� � ����������� ��������� ������� (� � ������� ����������� �� ����� �����,
			// ������� ����� ��������� � ������ �������.
			rLeonard_stress[i][j]=speed_test_filtering[fmr][j]*speed_test_filtering[smr][j]-speedxspeed_test_filtering[i][j];  
		}
	}
	// ������������ ������ �� ��� �������������
	// �������� ������������ ��������:
	if (speedxspeed_test_filtering != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (speedxspeed_test_filtering[i] != nullptr) {
				delete[] speedxspeed_test_filtering[i];
			}
		}
		delete[] speedxspeed_test_filtering;
	}

	if (speed_test_filtering != nullptr) {
		for (integer i = 0; i < 3; i++) {
			if (speed_test_filtering[i] != nullptr) {
				delete[] speed_test_filtering[i];
			}
		}
		delete[] speed_test_filtering;
	}

    // ����� ��������� ���������� ������� ��������� ���������� �� ������ 
	// ��� ����������� �������� ������� ����������� �� ��������� �������� �� ������� ����� ������.
	// ����� ���������� ������� ��������� ���������� ����� ���������.
	// ������ 9 �������� ��� ���������� ������� ��������� ����������.
	// ��������� �������� � �������� 9 ��� ������ ������� ��������� ����������.
	doublereal** StRT = nullptr;
	StRT = new doublereal*[10]; // Strain Rate Tensor
	for (integer i=0; i<10; i++) StRT[i]=new doublereal[f.maxelm+f.maxbound];
	const integer MODULEStRt=9; 

	// ���������� ��������� ������� �������� ����������.

	// ���� �� ���� ����������� �������
    #pragma omp parallel for
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		StRT[VXVX][i]=f.potent[GRADXVX][i];
		StRT[VXVY][i]=0.5*(f.potent[GRADYVX][i]+f.potent[GRADXVY][i]);
		StRT[VXVZ][i]=0.5*(f.potent[GRADZVX][i]+f.potent[GRADXVZ][i]);

		StRT[VYVX][i]=0.5*(f.potent[GRADXVY][i]+f.potent[GRADYVX][i]);
		StRT[VYVY][i]=f.potent[GRADYVY][i];
		StRT[VYVZ][i]=0.5*(f.potent[GRADZVY][i]+f.potent[GRADYVZ][i]);

		StRT[VZVX][i]=0.5*(f.potent[GRADXVZ][i]+f.potent[GRADZVX][i]);
		StRT[VZVY][i]=0.5*(f.potent[GRADYVZ][i]+f.potent[GRADZVY][i]);
		StRT[VZVZ][i]=f.potent[GRADZVZ][i];

		// ������ ������� ��������� ����������:
		StRT[MODULEStRt][i]=sqrt(2.0*(StRT[VXVX][i]*StRT[VXVX][i]+
			                     StRT[VYVY][i]*StRT[VYVY][i]+
								 StRT[VZVZ][i]*StRT[VZVZ][i])+
								 4.0*(StRT[VXVY][i]*StRT[VXVY][i]+
								 StRT[VXVZ][i]*StRT[VXVZ][i]+
								 StRT[VZVY][i]*StRT[VZVY][i]));
	}

	if (speedxspeed_buf != nullptr) {
		delete[] speedxspeed_buf;
		speedxspeed_buf = nullptr;
	}
	doublereal* Mij_buf = nullptr;
	Mij_buf = new doublereal[f.maxelm + f.maxbound];

	doublereal** Mij=new doublereal*[9]; // ������� ������������ Mij ������� ���������� ������.
	for (integer i=0; i<9; i++) Mij[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<9; i++) {
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			Mij[i][j]=0.0; // �������������
		}
	}


	for (integer i=0; i<9; i++) {

		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			Mij_buf[j]=StRT[MODULEStRt][j]*StRT[i][j];
		}

        // ���������� ������������� �������� ������������ ���������� � Mij_buf.
	    // �������� ! ������ �� ����������� ��������� �������� �� ����� ����������� �� ���� ������������� ��������.
	    // ��� ���������� �������� �� ������� ������� ����������� ������������ ������������ ������� ������� ������.
	    double_average_potent(Mij_buf, Mij[i], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	}

	if (Mij_buf != nullptr) {
		delete[] Mij_buf;
	}

	doublereal** StRT_filtr = nullptr;
	StRT_filtr = new doublereal*[10]; // ������ ������������� Strain Rate Tensor
	for (integer i=0; i<10; i++) StRT_filtr[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<10; i++) for (integer j=0; j<f.maxelm+f.maxbound; j++) StRT_filtr[i][j]=0.0; // �������������


    for (integer i=0; i<9; i++) {

		// ���������� ������ ��������������� ��������� Strain Rate Tensor`�.
	    // �������� ! ������ �� ����������� ��������� �������� �� ����� ����������� �� ���� ������������� ��������.
	    // ��� ���������� �������� �� ������� ������� ����������� ������������ ������������ ������� ������� ������.
	    double_average_potent(StRT[i], StRT_filtr[i], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	}

	delta_test_filtr=new doublereal[f.maxelm]; // ������ ��������� �������.
	for (integer i=0; i<f.maxelm; i++) delta_test_filtr[i]=0.0;

    // ���������� ������ ��������������� ��������� Strain Rate Tensor`�.
	// �������� ! ������ �� ����������� ��������� �������� �� ����� ����������� �� ���� ������������� ��������.
	// ��� ���������� �������� �� ������� ������� ����������� ������������ ������������ ������� ������� ������.
	double_average_potent(StRT[MODULEStRt], StRT_filtr[MODULEStRt], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	if (StRT != nullptr) {
		for (integer i = 0; i < 10; i++) {
			if (StRT[i] != nullptr) {
				delete[] StRT[i];
			}
		}
		delete[] StRT;
	}

	// �������-�� ���������� Mij:
	for (integer i=0; i<9; i++) {
		for (integer j=0; j<f.maxelm; j++) {
			doublereal dx=0.0, dy=0.0, dz=0.0; // ���������� ����� ��������
		    volume3D(j, f.nvtx, f.pa, dx, dy, dz);
		    doublereal rsvol=dx*dy*dz;
			doublereal delta2=exp((2.0/3.0)*log(rsvol));
			// ���� ������� ��� �� ��������������� ����������� ������ ������� �� ����� ���� ����������� �����.
			// ��� ����������� ����� �������� ������ ����� � ��� ���� ������ ��� ������� ������.
			// Mij[i][j]=(delta_test_filtr[j]*delta_test_filtr[j])*StRT_filtr[MODULEStRt][j]*StRT_filtr[i][j]-delta2*Mij[i][j];
			// ������� ������ ��������� ������� � 9 ��� ������ ������ �������� �������� �������.
			// ���� ������� � ������ ������ ����� ���� ����� �� ����� ���� ����������� �����.
			// �� ����������� ���������� ������������� �� ��������� ������� ����� ����� ����������.
			// ����� �� ����� ���� ����������� ����� ����������� �� ��������������� ����������� ������ 
			// ������������� ���������� �������������
            Mij[i][j]=delta2*(9.0*StRT_filtr[MODULEStRt][j]*StRT_filtr[i][j]-Mij[i][j]);
			// ���� ������ �� �������� �� �������� ������� � ��������� ����� ���������� ��������� ������������.
			// ���� ����� ���������� ��������� ������������ ������� � ��� ��� ����� �������� ����������� ������
			// (��������������� ������� ����������� �����) 26 ����������� � ��� ������������ ��������. ����� �������
			// ����� �������� 26 �������������� �������������� ������� � ������� ����� ������������ ����������� �������
			// �������������� ������� � ������������� ����� �� ����������� ����� ������� � ������� ������������ ������������.
			// ������ �� ��� �� ����� ����� ��� �� ��������� ������� ��������� ���� �������� ����������� ����� � �� ����������������
			// ����������� ������ ��������� �� ����� ����������� ����� ����������� �� �����.
			// ������: ��������� ���������� �������� ������������.
		}
	}

	// ������������ ����������� ������.

	if (StRT_filtr != nullptr) {
		for (integer i = 0; i < 10; i++) {
			if (StRT_filtr[i] != nullptr) {
				delete[] StRT_filtr[i];
			}
		}
		delete[] StRT_filtr;
	}

	if (delta_test_filtr != nullptr) {
		delete[] delta_test_filtr;
	}

	Cs2=new doublereal[f.maxelm];
	for (integer i=0; i<f.maxelm; i++) Cs2[i]=0.0; // �������������.

	// Cs2 == Cg - ��������� �������, ������� ��������� �������������.
	for (integer i=0; i<f.maxelm; i++) {
		doublereal rLijMij=0.0;
		for (integer j=0; j<9; j++) rLijMij+=rLeonard_stress[j][i]*Mij[j][i];
		doublereal rMijMij=0.0;
		for (integer j=0; j<9; j++) rMijMij+=Mij[j][i]*Mij[j][i];
		bool bzero=false; // ���� true �� ��������� ������� ����� ���� ��� �� ���� �������� ������ ����� �������.

		if (fabs(rMijMij)<1e-23) {
			if (fabs(rLijMij)<1e-23) {
				bzero=true;
			}
			else {
				printf("Error in Germano Dynamic Model: Mij is equal 0.0...\n");
			    printf("division by zero! ... \n");
			    printf("Please, press any key to exit...");
			    //getchar();
				system("pause");
			    exit(0);
			}
		}

		// ���������: �������� ��������� ������� ����� ������ ����������
		// � ������������ � �������, ��� ����� �������� � ��������� ��������������.
		// ������������ ���� ��������� � ��� ����������, ���� ��������� � ��� ������ 
		// ����������. 
		// � ��������� ������� ������������ ���������� ����� ��� ��������� �������������.
		// � [3] ������������ ����� 0.06 <= Cs <= 0.25. ��� ��� ��� �������� ��������� 
		// ������������� ��������: 0.0036 <= ��������� ������� <= 0.0625. ���������� 
		// ������� ������������ ������������ ��������������� ������ ��������������.
		// [3]. �.�. ���������, �.�. ��������, �.�.������ ��������� ������������� �������������
		// ������� � ������. ��������� ������������ ����������� ��. ���-������, ������.
		// e-mail: uali:kazsu.kz


		// �������� ����� ��� ����� �� ����� �����, 
		// ������ ����������� ���������� ��������.
		// ��� ��� ������� ������ ���� �����������.
		if (!bzero) {
			Cs2[i]=rLijMij/(2.0*rMijMij);

			// ������������ ����������� ��������������� ������.
		    //Cs2[i]=fmin(0.0036, fmax(Cs2[i], 0.0625)); 
		    if (f.smaginfo.bLimiters_Cs) {
				// ������������ ��������� �������������:
			    doublereal rfmax=f.smaginfo.maxCs*f.smaginfo.maxCs;
			    if (f.smaginfo.maxCs<0.0) rfmax*=-1.0;
			    doublereal rfmin=f.smaginfo.minCs*f.smaginfo.minCs;
			    if (f.smaginfo.minCs<0.0) rfmin*=-1.0;
			    if (rfmax>rfmin) {
                    //Cs2[i]=fmax( rfmin, fmin(Cs2[i], rfmax)); 
			    }
			    else {
				    printf("improper restriction of the constant Smagorinsky...\n");
				    printf("Please, press any key to exit is programm...\n");
				   // getchar();
					system("pause");
				    exit(0);
			    }
		    }
		    // ���� ������� ����������� ��� �� ������������ ����������� �������� �������� ���������� �������������.
		}
		else {
			// � ��������� � ����������� ����� ����, ������
			// ����� ������� ���� ��������, � ������ ���� �������
			// ��������������.
			Cs2[i]=0.0; 
		}
		
	}

	// � ���������� ����� ��������� �������� ���� ��������,
	// ��� ������������ �� ���������� ���������.
	// ������� ��������� ������� ������.
	// ����. ��������� � ��������� ��� ����������� �����.
	// � ����� ��� �� ������� � ������ a � ������ b.
	// ����� ��� ������������� b>>a. ����� ����� ��������� 
	// ������ �������� ������� � ������ ��������� ������� ����� 
	// ������� ��� ������ ���������������� ����������� �������.
	// ����� ��� ������ ������������� b==2.0*a;
	// �����: a/(3.0*a)==1.0/3.0, // ��� ����������� ����� ������ a.
	// a/(2.0*a+b)=a/(4.0*a)==1.0/4.0
	// b/(2.0*b+a)=2.0/(4.0+1.0)=2.0/5.0
	// b/(3.0*b)=1.0/3.0 // ��� ����������� �����.
	// �����:
	// a      a  < b   b 
	// 0.333 0.25 0.4 0.333
	// �� ����� ������� �� ����� �����.
	/*
	doublereal* Cs_flat=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		if (i<f.maxelm) {
	     	Cs_flat[i]=Cs2[i];
		}
		else Cs_flat[i]=0.0;
	}
	// �������� ������������ ���������� �������� �������������.
	// ������������ ����� ��������� �� ������������ �����������.
	flattener(Cs_flat, f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, f.nvtx, f.pa, f.border_neighbor);
	for (integer i=0; i<f.maxelm; i++) Cs2[i]=Cs_flat[i];
	delete Cs_flat;
	*/

	// ������������ ����������� ������.
	
	// ����������� ������ �� ��� ���������� ��������.
	if (rLeonard_stress != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (rLeonard_stress[i] != nullptr) {
				delete[] rLeonard_stress[i];
			}
		}
		delete[] rLeonard_stress;
	}

	if (imarker != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (imarker[i] != nullptr) {
				delete[] imarker[i];
			}
		}
		delete[] imarker;
	}

	if (Mij != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (Mij[i] != nullptr) {
				delete[] Mij[i];
			}
		}
		delete[] Mij;
	}

	

} // ������ �������.


#endif
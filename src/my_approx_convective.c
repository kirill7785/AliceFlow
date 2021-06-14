// ���� my_approx_convective.c 
// ������������� ������������� �����

#ifndef MY_APPROX_CONVECTIVE_C
#define MY_APPROX_CONVECTIVE_C 1

//#include <math.h>

//#define doublereal double

/*
*	�.�. ��������, �.�. �������� � ���������� ������ ������������ �����
* ������������� � �����������, ���������� �� ������ ������������ ������.
* ����������� ��������������� ����������� �����������. 2009, 12 �.
*/



// ���������� ������������ �� ���� 
// ������������ �����.
/*
doublereal fmax(doublereal fA, doublereal fB) {
	doublereal r=fB;
	if (fA > fB) r=fA;
	return r;
} // fmax
*/

// ������� A(|P|) 
// ��� ��������� ����:
// fabsPe - ������ ����� �����.
// ishconvection - ����� �����:
// CR==1 - ����������-���������� �����,
// UDS==2 - � ���������� ������ ������,
// COMB==3 - ���������������,
// POLY==4 - �� ��������� �������,
// EXP==5 - ���������������� (������),
// BULG==6 - ����� �.�. ��������� (23) �� ������,
// POW==7 - ������������� �����������.
doublereal ApproxConvective(doublereal fabsPe, integer ishconvection) {
	doublereal r=1; // �� ��������� � ��������� ������ ������
	doublereal rCR,  rpoly5;
	if (ishconvection<6) {
 	   rCR=1.0-0.5*fabsPe;
	   doublereal rpoly, rpoly2;
	   rpoly=1.0-0.1*fabsPe;
	   rpoly2=rpoly*rpoly;
	   rpoly5=rpoly2*rpoly2*rpoly; // �� ��������� �������
	}

	switch (ishconvection) {
		case CR: // ���������� ���������� ����� �������� ������ 
			    // ��� ��������� ����� ����� ������� �� ������ 1.0.
				if (fabsPe < 1.0) {
					r = rCR;
				}
				else {
					// ������� ��� ������ ����� �� ������ ������� 1.0 ��
					// ������� 10.0 �������������� ������� �� �������������
					// �. ���������� ����� �� ��������� �������.
					if ((fabsPe >= 1.0) && (fabsPe < 10.0)) r = rpoly5;
					else r = 0.0; // ��� ������ ����� ������� 10.0 ������������� ����.
					
				}
			    break;
		case UDS: // � ���������� ������ ������
			    //if (fabsPe<1.0) r=1.0; else r=0.0;
				r=1.0;
				//getchar();
			    break;
		case COMB: // ���������������
			    r=fmax(0.0, rCR);
			    break;
		case POLY: // �� ��������� �������
			    if (fabsPe<10.0) r=rpoly5; else r=0.0;
			    break;
		case EXP: // ���������������� (������)
			    if (fabsPe<10.0) {
			    	if (fabsPe<0.01) r=rCR; else r=fabsPe/expm1(fabsPe);
			    }
			    else r=0.0;
			    break;
		case BULG: // �.�. ��������, �.�. ��������
			    // ���������� ��������������� ����������� �����������
			    // ����������� (23) �� ������.
			    r=1.0/(1.0+0.6712*fabsPe*fabsPe);
			    break;
		case POW: // ������������� �����������
			    r=pow(0.553,fabsPe);
			    break;
		default: r=1.0/(1.0+0.6712*fabsPe*fabsPe); break;
	}
	return r;
} // ApproxConvective

#endif
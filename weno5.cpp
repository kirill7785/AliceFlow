

#pragma once
#ifndef MY_WENO5_CPP
#define MY_WENO5_CPP 1

typedef struct TValue_of_face {
	doublereal forvard, backward;
} Value_of_face;


// Схема WENO5 написана 13-14 января 2021
// для неравномерной сетки.
// Observations on the fifth-order WENO method with non-uniform meshes.
// Rong Wang, Hui Feng, Raymond J. Spiteri
Value_of_face weno5(
	doublereal h1,
	doublereal h2, 
	doublereal h3, 
	doublereal h4, 
	doublereal h5,
	doublereal u1,
	doublereal u2,
	doublereal u3,
	doublereal u4,
	doublereal u5)
{
	doublereal b22k = 1.0 / (h1 + h2 + h3) + 1 / (h2 + h3) + 1 / h3;
	doublereal b21k = b22k - ((h1 + h2 + h3) * (h2 + h3)) / ((h1+h2)*h2*h3);
	doublereal b20k = b21k + ((h1+h2+h3)*h3) / (h1*h2*(h2+h3));

	doublereal b12k = ((h2 + h3) * h3) / ((h2 + h3 + h4) * (h3 + h4) * h4);
	doublereal b11k = b12k + 1.0 / (h2 + h3) + 1.0 / h3 - 1.0 / h4;
	doublereal b10k = b11k - ((h2 + h3) * h4) / (h2*h3*(h3+h4));

	doublereal b02k = -h3 * h4 / ((h3+h4+h5)*(h4+h5)*h5);
	doublereal b01k = b02k + (h3 * (h4 + h5)) / ((h3 + h4) * h4 * h5);
	doublereal b00k = b01k + 1.0 / h3 - 1.0 / h4 - 1.0 / (h4 + h5);

	doublereal b22 = b12k, b21 = b11k, b20 = b10k;
	doublereal b12 = b02k, b11 = b01k, b10 = b00k;
	doublereal b02 = (h3 * (h3 + h4)) / ((h3+h4+h5)*(h4+h5)*h5);
	doublereal b01 = b02 - (h3 * (h3 + h4 + h5)) / ((h3+h4)*h4*h5);
	doublereal b00 = b01 + ((h3+h4)*(h3+h4+h5)) / (h3*h4*(h4+h5));

	doublereal P2forvard = b20k * h1 * u1 + b21k * h2 * u2 + b22k * h3 * u3;
	doublereal P1forvard = b10k * h2 * u2 + b11k * h3 * u3 + b12k * h4 * u4;
	doublereal P0forvard = b00k * h3 * u3 + b01k * h4 * u4 + b02k * h5 * u5;

	doublereal P2backward = b20 * h1 * u1 + b21 * h2 * u2 + b22 * h3 * u3;
	doublereal P1backward = b10 * h2 * u2 + b11 * h3 * u3 + b12 * h4 * u4;
	doublereal P0backward = b00 * h3 * u3 + b01 * h4 * u4 + b02 * h5 * u5;

	doublereal d2 = ((h3 + h4) * (h3 + h4 + h5)) / ((h1+h2+h3+h4)*(h1 + h2 + h3 + h4+h5));
	doublereal d1 = ((h1+h2)*(h3+h4+h5)*(h1+2.0*h2+2.0*h3+2.0*h4+h5)) / ((h1 + h2 + h3 + h4)* (h2 + h3 + h4 + h5) * (h1 + h2 + h3 + h4 + h5));
	doublereal d0 = (h2 * (h1 + h2)) / ((h2 + h3 + h4 + h5)*(h1 + h2 + h3 + h4 + h5));

	doublereal d2k = (h4*(h4+h5))/ ((h1 + h2 + h3 + h4) * (h1 + h2 + h3 + h4 + h5));
	doublereal d1k = ((h1+h2+h3)*(h4+h5)*(h1 + 2.0 * h2 + 2.0 * h3 + 2.0 * h4 + h5))/ ((h1 + h2 + h3 + h4) * (h2 + h3 + h4 + h5) * (h1 + h2 + h3 + h4 + h5));
	doublereal d0k = ((h2+h3)*(h1+h2+h3))/((h2 + h3 + h4 + h5) * (h1 + h2 + h3 + h4 + h5));


	doublereal B02ss = 6 / ((h3 + h4 + h5) * (h4 + h5) * h5);
	doublereal B12ss = 6 / ((h2 + h3 + h4) * (h3 + h4) * h4);
	doublereal B22ss = 6 / ((h1 + h2 + h3) * (h2 + h3) * h3);

	doublereal B01ss = B02ss - 6 / ((h3 + h4) * h4 * h5);
	doublereal B11ss = B12ss - 6 / ((h2 + h3) * h3 * h4);
	doublereal B21ss = B22ss - 6 / ((h1 + h2) * h2 * h3);

	doublereal B00ss = B01ss + 6 / (h3 * h4 * (h4 + h5));
	doublereal B10ss = B11ss + 6 / (h2 * h3 * (h3 + h4));
	doublereal B20ss = B21ss + 6 / (h1 * h2 * (h2 + h3));

	doublereal B22sm = (2.0*(h1+2.0*h2)) / ((h1+h2+h3)*(h2+h3)*h3);
	doublereal B21sm = B22sm - (2*(h1+2*h2-h3)) / ((h1+h2)*h2*h3);
	doublereal B20sm = B21sm + (2*(h1+h2-h3)) / (h1*h2*(h2+h3));

	doublereal B12sm = (2.0 * (h2-h3)) / ((h2 + h3 + h4) * (h3 + h4) * h4);
	doublereal B11sm = B12sm - (2 * (h2 - h3 - h4)) / ((h2 + h3) * h3 * h4);
	doublereal B10sm = B11sm + (2 * (h2 -2*h3 - h4)) / (h2 * h3 * (h3 + h4));

	doublereal B02sm = (4.0 * (h3 +  h4)) / ((h3 + h4 + h5) * (h4 + h5) * h5);
	doublereal B01sm = B02sm + (2 * (2*h3 + h4 + h5)) / ((h3 + h4) * h4 * h5);
	doublereal B00sm = B01sm - (2 * (2*h3 + 2*h4 + h5)) / (h3 * h4 * (h4 + h5));

	doublereal B22sc = B22sm + 0.5 * h3 * B22ss;
	doublereal B12sc = B12sm + 0.5 * h3 * B12ss;
	doublereal B02sc = B02sm + 0.5 * h3 * B02ss;

	doublereal B21sc = B21sm + 0.5 * h3 * B21ss;
	doublereal B11sc = B11sm + 0.5 * h3 * B11ss;
	doublereal B01sc = B01sm + 0.5 * h3 * B01ss;

	doublereal B20sc = B20sm + 0.5 * h3 * B20ss;
	doublereal B10sc = B10sm + 0.5 * h3 * B10ss;
	doublereal B00sc = B00sm + 0.5 * h3 * B00ss;

	doublereal B22sp = B22sm + h3 * B22ss;
	doublereal B12sp = B12sm + h3 * B12ss;
	doublereal B02sp = B02sm + h3 * B02ss;

	doublereal B21sp = B21sm + h3 * B21ss;
	doublereal B11sp = B11sm + h3 * B11ss;
	doublereal B01sp = B01sm + h3 * B01ss;

	doublereal B20sp = B20sm + h3 * B20ss;
	doublereal B10sp = B10sm + h3 * B10ss;
	doublereal B00sp = B00sm + h3 * B00ss;

	doublereal IS2 = 0.0;
	doublereal IS1 = 0.0;
	doublereal IS0 = 0.0;


	IS2 = h3 * h3 * (pow((1.0/6.0)*(B20sm*h1*u1+ B21sm * h2 * u2 + B22sm * h3 * u3 ),2.0)+
		4.0*pow((1.0 / 6.0) * (B20sc*h1*u1+ B21sc * h2 * u2 + B22sc * h3 * u3 ),2.0)+
		pow((1.0 / 6.0) * (B20sp*h1*u1+ B21sp * h2 * u2 + B22sp * h3 * u3 ),2.0));

	IS1 = h3 * h3 * (pow((1.0 / 6.0) * (B10sm * h2 * u2 + B11sm * h3 * u3 + B12sm * h4 * u4), 2.0) +
		4.0 * pow((1.0 / 6.0) * (B10sc * h2 * u2 + B11sc * h3 * u3 + B12sc * h4 * u4), 2.0) +
		pow((1.0 / 6.0) * (B10sp * h2 * u2 + B11sp * h3 * u3 + B12sp * h4 * u4), 2.0));

	IS0 = h3 * h3 * (pow((1.0 / 6.0) * (B00sm * h3 * u3 + B01sm * h4 * u4 + B02sm * h5 * u5), 2.0) +
		4.0 * pow((1.0 / 6.0) * (B00sc * h3 * u3 + B01sc * h4 * u4 + B02sc * h5 * u5), 2.0) +
		pow((1.0 / 6.0) * (B00sp * h3 * u3 + B01sp * h4 * u4 + B02sp * h5 * u5), 2.0));

	IS2 += pow(h3, 4.0) * pow(B20ss * h1 * u1 + B21ss * h2 * u2 + B22ss * h3 * u3 ,2.0);
	IS1 += pow(h3, 4.0) * pow(B10ss * h2 * u2 + B11ss * h3 * u3 + B12ss * h4 * u4, 2.0);
	IS0 += pow(h3, 4.0) * pow(B00ss * h3 * u3 + B01ss * h4 * u4 + B02ss * h5 * u5, 2.0);

	const doublereal epsilon = 1.0e-6;

	doublereal alpha2k = d2k / ((epsilon + IS2) * (epsilon + IS2));
	doublereal alpha1k = d1k / ((epsilon + IS1) * (epsilon + IS1));
	doublereal alpha0k = d0k / ((epsilon + IS0) * (epsilon + IS0));

	doublereal alpha2 = d2 / ((epsilon + IS2) * (epsilon + IS2));
	doublereal alpha1 = d1 / ((epsilon + IS1) * (epsilon + IS1));
	doublereal alpha0 = d0 / ((epsilon + IS0) * (epsilon + IS0));

	doublereal w2k = alpha2k / (alpha0k + alpha1k + alpha2k);
	doublereal w1k = alpha1k / (alpha0k + alpha1k + alpha2k);
	doublereal w0k = alpha0k / (alpha0k + alpha1k + alpha2k);

	doublereal w2 = alpha2 / (alpha0 + alpha1 + alpha2);
	doublereal w1 = alpha1 / (alpha0 + alpha1 + alpha2);
	doublereal w0 = alpha0 / (alpha0 + alpha1 + alpha2);

	Value_of_face ff;

    ff.forvard = w2k * P2forvard + w1k * P1forvard + w0k * P0forvard;
	ff.backward = w2 * P2backward + w1 * P1backward + w0 * P0backward;

	return ff;

}// weno5

#endif
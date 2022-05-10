// ChebyshevSmoother.cpp : 
// Пафнутий Львович Чебышев 4 (16) мая 1821 - 26 ноября (8 декабря) 1894.
// 6.01.2022 Пафнутию Львовичу Чебышеву 200 лет.
// В данном файле содержится реализация сглаживателя Пафнутия Львовича Чебышева для алгебраического многосеточного метода РУМБА v.0.14
// и сглаживатель Чебышева применим только для положительно определенных симметричных задач (Механика, теплопередача в твёрдом теле, поправка давления для гидродинамики.
// Особенностью данного файла является следующее:
// 1. работает для трех типов хранения разреженной матрицы 1.1. CRS (val, col_ind, row_ptr) делается выделение и уничножение памяти, используется дополнительная память под еще одну данную матрицу.
// 1.2. AK1* &Amat (не протестировано), 1.3. Ak2 &Amat (протестировано). Причём Ak1 и Ak2 работают на месте и лишнего выделения памяти под новую матрицу не делаю. 
// Работа на месте даёт ускорение в время расчёта 8мин 27с 470мс против 10мин 21с 320мс с выделениями и освобождеениями памяти. На примере Механики для фермы, 137тысяч неизвестных, степень полинома Чебышева 64.
// Также здесь в файле ChebyshevSmoother.cpp присутствует вариант сглаживателя Чебышева Дениса Демидова по мотивам библиотеки AMGCL. Заимствованый у ddemidova вариант работает только с матрицей  CRS (val, col_ind, row_ptr).
// Кроме варианта ddemidovа присутствует адаптивный двухслойный вариант сглаживателя Чебышева у которого нижняя граница спектра уточняется - материалам Жукова, Феодоритова для адаптации нижней границы спектра, и сами 
// оптмальные Чебышевские параметры для двухслойного итерационного метода по мтивам книги Самарского и Николаева. 
// Используется переупорядочивание корней многочлена Чебышева для устойчивости вычислительного алгоритма (впервые Янг 1954).
// Точное решение задачи переупорядочивания нулей многочлена Чебышева В.И. Лебедев 1970 годы.

// Ричардсоном в работе, опубликованной в 1910 году, был предложен способ ускорения сходимости. В этом способе итерационный процесс имел также свои параметры на каждом шаге,
// но выбор параметров не был оптимальным, хотя и был связан со свойствами некоторых алгебрических многочленов. Требования оптимальности обеспечивают полиномы Чебышева,
// но Ричардсон их не привел в своей работе, несмотря на широкую известность полиномов Чебышева к тому времени.
// Оптимальные многочлены были построены П.Л.Чебышевым уже в середине XIX века, затем исследования были продолжены его учениками Е.И. Золотарёвым и А.А.Марковым. 

#pragma once
#ifndef BASIC_FUNCTIONS_CHEBYSHEV_SMOOTHERS_CPP
#define BASIC_FUNCTIONS__CHEBYSHEV_SMOOTHERS_CPP 1

#include <iostream>
#include <random>

#include "cache.cpp"


// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.

doublereal spectral_radius(doublereal* val, integer* col_ind, integer* row_ptr, integer n, bool scale, integer power_iters = 0) {
    //std::cout << "spectral radius" << std::endl;

   
    doublereal radius;

    if (power_iters <= 0) {
        // Use Gershgorin disk theorem.
        radius = 0;

#pragma omp parallel
        {
            doublereal emax = 0;
            doublereal  dia = 1.0;

#pragma omp for nowait
            for (integer i = 0; i < n; ++i) {
                doublereal s = 0;

                for (integer j = row_ptr[i], e = row_ptr[i + 1]; j < e; ++j) {
                    integer  c = col_ind[j];
                    doublereal v = val[j];

                    s += fabs(v);// норма числа

                    if (scale && c == i) dia = v;
                }

                if (scale) s *= 1.0/(fabs(dia));

                emax = std::max(emax, s);
            }

#pragma omp critical
            radius = std::max(radius, emax);
        }
    }
    else {
        // Power method.
        
        doublereal* b0 = new doublereal[n];
        doublereal* b1 = new doublereal[n];
        //for (integer i = 0; i < n; ++i) {
           // b0[i] = 0.0;
           // b1[i] = 0.0;
        //}

        // Fill the initial vector with random values.
        // Also extract the inverted matrix diagonal values.
        doublereal b0_norm = 0.0;
#pragma omp parallel
        {
#ifdef _OPENMP
            unsigned int tid = omp_get_thread_num();
#else
            unsigned int tid = 0;
#endif
            std::mt19937 rng(tid);
            std::uniform_real_distribution<doublereal> rnd(-1, 1);

            doublereal loc_norm = 0.0;

#pragma omp for nowait
            for (integer i = 0; i < n; ++i) {
                doublereal v = (rnd(rng));

                b0[i] = v;
                loc_norm += v * v; // math::norm(math::inner_product(v, v));
            }

#pragma omp critical
            b0_norm += loc_norm;
        }

        // Normalize b0
        b0_norm = 1 / sqrt(b0_norm);
#pragma omp parallel for
        for (integer i = 0; i < n; ++i) {
            b0[i] = b0_norm * b0[i];
        }

        for (integer iter = 0; iter < power_iters;) {
            // b1 = scale ? (D^1 * A) * b0 : A * b0
            // b1_norm = ||b1||
            // radius = <b1,b0>
            doublereal b1_norm = 0.0;
            radius = 0.0;
#pragma omp parallel
            {
                doublereal loc_norm = 0.0;
                doublereal loc_radi = 0.0;
                doublereal  dia = 1.0;

#pragma omp for nowait
                for (integer i = 0; i < n; ++i) {
                    doublereal s = 0.0;

                    for (integer j = row_ptr[i], e = row_ptr[i + 1]; j < e; ++j) {
                        integer  c = col_ind[j];
                        doublereal v = val[j]; //val[j]; // моя фантазия fabs(val[j]); тоже не работает.
                        if (scale && c == i) dia = v;
                        s += v * b0[c];
                    }

                    if (scale) s = (1.0/(dia)) * s;

                    loc_norm += s * s;
                    loc_radi += fabs(s * b0[i]);

                    b1[i] = s;
                }

#pragma omp critical
                {
                    b1_norm += loc_norm;
                    radius += loc_radi;
                }
            }

            if (++iter < power_iters) {
                // b0 = b1 / b1_norm
                b1_norm = 1 / sqrt(b1_norm);
#pragma omp parallel for
                for (integer i = 0; i < n; ++i) {
                    b0[i] = b1_norm * b1[i];
                }
            }
        }
    }
    //std::cout << "spectral radius" << std::endl;

    return radius < 0 ? static_cast<doublereal>(2.0) : radius;
}

// Вычисление невязки (специально для метода Чебышева).
inline void residual_for_Cheb(doublereal* &val, integer* &col_ind, integer* &row_ptr, doublereal* &b, doublereal* &x, doublereal*& r, integer n)
{
    //MatrixCRSByVector(val, col_ind, row_ptr, x, r, n);

#pragma omp parallel for 
    for (integer i = 0; i < n; ++i) {
        doublereal sum = 0.0;
        const integer rowend = row_ptr[i + 1];
        const integer rowbeg = row_ptr[i];


        for (integer j = rowbeg; j < rowend; ++j)
        {
            sum += val[j] * x[col_ind[j]];
        }
        r[i] = sum;
    }

#pragma omp parallel for
    for (integer i = 0; i < n; ++i) r[i] = b[i] - r[i];

} // residual_for_Cheb

void chebyshev0( doublereal *val, integer *col_ind, integer *row_ptr, doublereal *b, doublereal *x, integer n, integer iVar, integer level)
{
    
    /// Chebyshev polynomial degree.
    integer degree;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    doublereal higher;

    /// Lowest-to-highest eigen value ratio.
    doublereal lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    integer power_iters;

    // Scale the system matrix
    bool scale;

    degree = my_amg_manager.Chebyshev_degree;
    higher = 1.0f;
    lower = 1.0f / 30;
    power_iters = 0;
    scale = false;
    
    
    doublereal hi, lo;


    // Спектральный радиус вычисляем если он еще не был вычислен ранее по теореме Гершгорина.
    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Temp) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Temp = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Temp = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Temp = true;
        }
        hi = level_Chebyshev_info[level].hi_Temp;
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Pressure) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Pressure = true;
        }
        hi = level_Chebyshev_info[level].hi_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Stress) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Stress = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Stress = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Stress = true;
        }
        hi = level_Chebyshev_info[level].hi_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vx) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vx = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vx = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vx = true;
        }
        hi = level_Chebyshev_info[level].hi_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vy) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vy = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vy = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vy = true;
        }
        hi = level_Chebyshev_info[level].hi_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vz) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vz = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vz = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vz = true;
        }
        hi = level_Chebyshev_info[level].hi_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_nu) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_nu = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_nu = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_nu = true;
        }
        hi = level_Chebyshev_info[level].hi_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_k = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k = true;
        }
        hi = level_Chebyshev_info[level].hi_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_omega) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_omega = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_omega = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_omega = true;
        }
        hi = level_Chebyshev_info[level].hi_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k_for_ke) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k_for_ke = true;
        }
        hi = level_Chebyshev_info[level].hi_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_epsilon) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_epsilon = true;
        }
        hi = level_Chebyshev_info[level].hi_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_ReTheta) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_ReTheta = true;
        }
        hi = level_Chebyshev_info[level].hi_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_gamma) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_gamma = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_gamma = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_gamma = true;
        }
        hi = level_Chebyshev_info[level].hi_gamma;
        break;
    default:

        if (!level_Chebyshev_info[level].bGershgorin_calc_Speed) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Speed = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Speed = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Speed = true;
        }
        hi = level_Chebyshev_info[level].hi_Speed;
        break;
    }


    lo = hi * lower;
    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();

    doublereal alpha = 0.0, beta = 0.0;

    doublereal* r = new doublereal[n];
    doublereal* p = new doublereal[n];

    // Инициализация нулевым значением.
#pragma omp parallel for
    for (integer i = 0; i < n; ++i) {
        r[i] = 0.0;
        p[i] = 0.0;
    }

    for (integer k = 0; k < degree; ++k) {
        //backend::residual(b, A, x, *r);
        residual_for_Cheb(val, col_ind, row_ptr, b, x, r, n);

        //if (prm.scale) backend::vmul(one, *M, *r, zero, *r);

        if (k == 0) {
            alpha = 1.0/d;
            beta = 0.0;
        }
        else if (k == 1) {
            alpha = 2 * d * (1.0/(2 * d * d - c * c));
            beta = alpha * d - 1.0;
        }
        else {
            alpha = 1.0 / ((d - 0.25 * alpha * c * c));
            beta = alpha * d - 1.0;
        }

        if (fabs(beta) > 1.0e-30) {
#pragma omp parallel for
            for (integer i = 0; i < n; ++i) {
                p[i] = alpha * r[i] + beta * p[i];
                x[i] += p[i];
            }
        }
        else {
#pragma omp parallel for
            for (integer i = 0; i < n; ++i) {
                p[i] = alpha * r[i];
                x[i] += p[i];
            }
        }
        
    }

    delete[] r;
    delete[] p;
}




// Уточнение нижней границы спектра.
void chebyshev_ddemidov(doublereal *val, integer *col_ind, integer *row_ptr, doublereal *b, doublereal *x, integer n, integer iVar, integer level)
{

    /// Chebyshev polynomial degree.
    integer degree;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    integer power_iters;

    // Scale the system matrix
    bool scale;

    degree = my_amg_manager.Chebyshev_degree;
    higher = 1.0f;
    lower = 1.0f / 30;
    power_iters = 0;
    scale = false;


    doublereal hi, lo;

    

    // Спектральный радиус вычисляем если он еще не был вычислен ранее по теореме Гершгорина.
    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Temp) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Temp = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Temp = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Temp = true;
        }
        hi = level_Chebyshev_info[level].hi_Temp;
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Pressure) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Pressure = true;
        }
        hi = level_Chebyshev_info[level].hi_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Stress) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Stress = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Stress = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Stress = true;
        }
        hi = level_Chebyshev_info[level].hi_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vx) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vx = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vx = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vx = true;
        }
        hi = level_Chebyshev_info[level].hi_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vy) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vy = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vy = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vy = true;
        }
        hi = level_Chebyshev_info[level].hi_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vz) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vz = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vz = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vz = true;
        }
        hi = level_Chebyshev_info[level].hi_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_nu) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_nu = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_nu = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_nu = true;
        }
        hi = level_Chebyshev_info[level].hi_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_k = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k = true;
        }
        hi = level_Chebyshev_info[level].hi_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_omega) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_omega = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_omega = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_omega = true;
        }
        hi = level_Chebyshev_info[level].hi_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k_for_ke) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k_for_ke = true;
        }
        hi = level_Chebyshev_info[level].hi_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_epsilon) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_epsilon = true;
        }
        hi = level_Chebyshev_info[level].hi_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_ReTheta) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_ReTheta = true;
        }
        hi = level_Chebyshev_info[level].hi_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_gamma) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_gamma = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_gamma = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_gamma = true;
        }
        hi = level_Chebyshev_info[level].hi_gamma;
        break;
    default:

        if (!level_Chebyshev_info[level].bGershgorin_calc_Speed) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Speed = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Speed = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Speed = true;
        }
        hi = level_Chebyshev_info[level].hi_Speed;
        break;
    }



    doublereal* r = new doublereal[n];
    doublereal* p = new doublereal[n];

    // Инициализация нулевым значением.
#pragma omp parallel for
    for (integer i = 0; i < n; ++i) {
        r[i] = 0.0;
        p[i] = 0.0;
    }


    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bstart_lo_Temp) {
            level_Chebyshev_info[level].lo_Temp  = hi * lower;

    // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Temp = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Temp = true;
        lo = level_Chebyshev_info[level].lo_Temp;
        //if (lo < 1.0e-14)
        {
            //std::cout << "incomming lo = " << lo << " level=" << level << std::endl;
          //  getchar();
        }
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bstart_lo_Pressure) {
            level_Chebyshev_info[level].lo_Pressure  = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Pressure = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Pressure = true;
        lo = level_Chebyshev_info[level].lo_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bstart_lo_Stress) {
            level_Chebyshev_info[level].lo_Stress  = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Stress = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Stress = true;
        lo = level_Chebyshev_info[level].lo_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vx) {
            level_Chebyshev_info[level].lo_Vx  = hi * lower;

    // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Vx = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vx = true;
        lo = level_Chebyshev_info[level].lo_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vy) {
            level_Chebyshev_info[level].lo_Vy  = hi * lower;

    // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Vy = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vy = true;
        lo = level_Chebyshev_info[level].lo_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vz) {
            level_Chebyshev_info[level].lo_Vz  = hi * lower;

    // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Vz = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vz = true;
        lo = level_Chebyshev_info[level].lo_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bstart_lo_nu) {
            level_Chebyshev_info[level].lo_nu = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_nu = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_nu = true;
        lo = level_Chebyshev_info[level].lo_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bstart_lo_k) {
            level_Chebyshev_info[level].lo_k = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k = true;
        lo = level_Chebyshev_info[level].lo_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bstart_lo_omega) {
            level_Chebyshev_info[level].lo_omega = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_omega = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_omega = true;
        lo = level_Chebyshev_info[level].lo_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_k_for_ke) {
            level_Chebyshev_info[level].lo_k_for_ke = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k_for_ke = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k_for_ke = true;
        lo = level_Chebyshev_info[level].lo_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_epsilon) {
            level_Chebyshev_info[level].lo_epsilon = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_epsilon = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_epsilon = true;
        lo = level_Chebyshev_info[level].lo_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_ReTheta) {
            level_Chebyshev_info[level].lo_ReTheta = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_ReTheta = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_ReTheta = true;
        lo = level_Chebyshev_info[level].lo_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_gamma) {
            level_Chebyshev_info[level].lo_gamma = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_gamma = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_gamma = true;
        lo = level_Chebyshev_info[level].lo_gamma;
        break;
    default:
        if (!level_Chebyshev_info[level].bstart_lo_Speed) {
            level_Chebyshev_info[level].lo_Speed  = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Speed = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Speed = true;
        lo = level_Chebyshev_info[level].lo_Speed;
        break;
    }



    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();

    doublereal alpha = 0.0, beta = 0.0;

    doublereal r_start=1.0;

    for (unsigned k = 0; k < degree; ++k) {
        //backend::residual(b, A, x, *r);
        residual_for_Cheb(val, col_ind, row_ptr, b, x, r, n);

        if (k == 0) r_start = NormaV(r, n);

        //if (prm.scale) backend::vmul(one, *M, *r, zero, *r);

        if (k == 0) {
            alpha = 1.0 / d;
            beta = 0.0;
        }
        else if (k == 1) {
            alpha = 2 * d * (1.0 / (2 * d * d - c * c));
            beta = alpha * d - 1.0;
        }
        else {
            alpha = 1.0 / ((d - 0.25 * alpha * c * c));
            beta = alpha * d - 1.0;
        }

        if (fabs(beta) > 1.0e-30) {
#pragma omp parallel for
            for (integer i = 0; i < n; ++i) {
                p[i] = alpha * r[i] + beta * p[i];
                x[i] += p[i];
            }
        }
        else {
#pragma omp parallel for
            for (integer i = 0; i < n; ++i) {
                p[i] = alpha * r[i];
                x[i] += p[i];
            }
        }

    }

    // calculate delta
    residual_for_Cheb(val, col_ind, row_ptr, b, x, r, n);

    doublereal delta = NormaV(r, n) / r_start;

    delete[] r;
    delete[] p;


    
   

    
    // Итерационное уточнение нижней границы спектра.
    // 1.
    doublereal etta = lo / hi;
    doublereal rho1 = (1.0 - sqrt(etta)) / (1.0 + sqrt(etta));
    doublereal qp = 2.0 * pow(rho1, 1.0*degree) / (1.0 + pow(rho1, 2.0 * degree));
    // 2.
    doublereal y1 = delta / qp;
    doublereal y2 = log(y1 + sqrt(  y1 * y1 - 1.0));
    // 3.
    doublereal x_zv = cosh(y2 / (1.0 * degree));
    // 4.
    doublereal etta_new =  0.5 * (1.0 + etta) - 0.5 * (1.0 - etta) * x_zv;
    // 5. 
    //std::cout << "x_zv=" << x_zv << " y2=" << y2 << " y1=" << y1 << " qp=" << qp << " rho1=" << rho1 << " etta=" << etta << " delta=" << delta << std::endl;
    //std::cout << "apostoriory lo = " << lo << "  etta_new=" << etta_new << "  level=" << level << std::endl;
    // getchar();
    

    if (etta_new > 1.0e-20) {

        switch (iVar) {
        case TEMP:
            level_Chebyshev_info[level].lo_Temp = etta_new * level_Chebyshev_info[level].hi_Temp;
            break;
        case PAM:
            level_Chebyshev_info[level].lo_Pressure = etta_new * level_Chebyshev_info[level].hi_Pressure;
            break;
        case TOTALDEFORMATIONVAR:
            level_Chebyshev_info[level].lo_Stress = etta_new * level_Chebyshev_info[level].hi_Stress;
            break;
        case VELOCITY_X_COMPONENT:
            level_Chebyshev_info[level].lo_Vx = etta_new * level_Chebyshev_info[level].hi_Vx;
            break;
        case VELOCITY_Y_COMPONENT:
            level_Chebyshev_info[level].lo_Vy = etta_new * level_Chebyshev_info[level].hi_Vy;
            break;
        case VELOCITY_Z_COMPONENT:
            level_Chebyshev_info[level].lo_Vz = etta_new * level_Chebyshev_info[level].hi_Vz;
            break;
        case NUSHA:
            level_Chebyshev_info[level].lo_nu = etta_new * level_Chebyshev_info[level].hi_nu;
            break;
        case TURBULENT_KINETIK_ENERGY:
            level_Chebyshev_info[level].lo_k = etta_new * level_Chebyshev_info[level].hi_k;
            break;
        case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
            level_Chebyshev_info[level].lo_omega = etta_new * level_Chebyshev_info[level].hi_omega;
            break;
        case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
            level_Chebyshev_info[level].lo_k_for_ke = etta_new * level_Chebyshev_info[level].hi_k_for_ke;
            break;
        case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
            level_Chebyshev_info[level].lo_epsilon = etta_new * level_Chebyshev_info[level].hi_epsilon;
            break;
        case RE_THETA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_ReTheta = etta_new * level_Chebyshev_info[level].hi_ReTheta;
            break;
        case GAMMA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_gamma = etta_new * level_Chebyshev_info[level].hi_gamma;
            break;
        default:
            level_Chebyshev_info[level].lo_Speed = etta_new * level_Chebyshev_info[level].hi_Speed;
            break;
        }
    }
    
}


// Уточнение нижней границы спектра.
void chebyshev(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* b, doublereal* x, integer n, integer iVar, integer level)
{
    // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.


    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    integer power_iters;

    // Scale the system matrix
    bool scale;

    
    higher = 1.0f;
    lower = 1.0f / 30;
    power_iters = 0;
    scale = false;


    doublereal hi, lo;



    // Спектральный радиус вычисляем если он еще не был вычислен ранее по теореме Гершгорина.
    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Temp) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Temp = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Temp = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Temp = true;
        }
        hi = level_Chebyshev_info[level].hi_Temp;
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Pressure) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Pressure = true;
        }
        hi = level_Chebyshev_info[level].hi_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Stress) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Stress = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Stress = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Stress = true;
        }
        hi = level_Chebyshev_info[level].hi_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vx) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vx = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vx = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vx = true;
        }
        hi = level_Chebyshev_info[level].hi_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vy) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vy = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vy = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vy = true;
        }
        hi = level_Chebyshev_info[level].hi_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vz) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vz = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Vz = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vz = true;
        }
        hi = level_Chebyshev_info[level].hi_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_nu) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_nu = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_nu = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_nu = true;
        }
        hi = level_Chebyshev_info[level].hi_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_k = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k = true;
        }
        hi = level_Chebyshev_info[level].hi_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_omega) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_omega = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_omega = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_omega = true;
        }
        hi = level_Chebyshev_info[level].hi_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k_for_ke) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k_for_ke = true;
        }
        hi = level_Chebyshev_info[level].hi_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_epsilon) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_epsilon = true;
        }
        hi = level_Chebyshev_info[level].hi_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_ReTheta) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_ReTheta = true;
        }
        hi = level_Chebyshev_info[level].hi_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_gamma) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_gamma = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_gamma = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_gamma = true;
        }
        hi = level_Chebyshev_info[level].hi_gamma;
        break;
    default:

        if (!level_Chebyshev_info[level].bGershgorin_calc_Speed) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Speed = spectral_radius(val, col_ind, row_ptr, n, true, power_iters);
            }
            else {
                level_Chebyshev_info[level].hi_Speed = spectral_radius(val, col_ind, row_ptr, n, false, power_iters);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Speed = true;
        }
        hi = level_Chebyshev_info[level].hi_Speed;
        break;
    }



    doublereal* r = new doublereal[n];
    


    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bstart_lo_Temp) {
            level_Chebyshev_info[level].lo_Temp = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Temp = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Temp = true;
        lo = level_Chebyshev_info[level].lo_Temp;
        //if (lo < 1.0e-14)
        {
            //std::cout << "incomming lo = " << lo << " level=" << level << std::endl;
          //  getchar();
        }
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bstart_lo_Pressure) {
            level_Chebyshev_info[level].lo_Pressure = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Pressure = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Pressure = true;
        lo = level_Chebyshev_info[level].lo_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bstart_lo_Stress) {
            level_Chebyshev_info[level].lo_Stress = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Stress = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Stress = true;
        lo = level_Chebyshev_info[level].lo_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vx) {
            level_Chebyshev_info[level].lo_Vx = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vx = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vx = true;
        lo = level_Chebyshev_info[level].lo_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vy) {
            level_Chebyshev_info[level].lo_Vy = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vy = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vy = true;
        lo = level_Chebyshev_info[level].lo_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vz) {
            level_Chebyshev_info[level].lo_Vz = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vz = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vz = true;
        lo = level_Chebyshev_info[level].lo_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bstart_lo_nu) {
            level_Chebyshev_info[level].lo_nu = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_nu = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_nu = true;
        lo = level_Chebyshev_info[level].lo_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bstart_lo_k) {
            level_Chebyshev_info[level].lo_k = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k = true;
        lo = level_Chebyshev_info[level].lo_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bstart_lo_omega) {
            level_Chebyshev_info[level].lo_omega = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_omega = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_omega = true;
        lo = level_Chebyshev_info[level].lo_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_k_for_ke) {
            level_Chebyshev_info[level].lo_k_for_ke = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k_for_ke = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k_for_ke = true;
        lo = level_Chebyshev_info[level].lo_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_epsilon) {
            level_Chebyshev_info[level].lo_epsilon = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_epsilon = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_epsilon = true;
        lo = level_Chebyshev_info[level].lo_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_ReTheta) {
            level_Chebyshev_info[level].lo_ReTheta = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_ReTheta = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_ReTheta = true;
        lo = level_Chebyshev_info[level].lo_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_gamma) {
            level_Chebyshev_info[level].lo_gamma = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_gamma = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_gamma = true;
        lo = level_Chebyshev_info[level].lo_gamma;
        break;
    default:
        if (!level_Chebyshev_info[level].bstart_lo_Speed) {
            level_Chebyshev_info[level].lo_Speed = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Speed = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Speed = true;
        lo = level_Chebyshev_info[level].lo_Speed;
        break;
    }



    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();

    

    doublereal r_start;

    const doublereal tau0 = 2.0 / (hi+lo);
    const doublereal ksi = lo / hi;
    const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

   
    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu=new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);

    for (unsigned k = 0; k < degree; ++k) {
        //backend::residual(b, A, x, *r);
        residual_for_Cheb(val, col_ind, row_ptr, b, x, r, n);

        if (k == 0) r_start = NormaV(r, n);

        //if (prm.scale) backend::vmul(one, *M, *r, zero, *r);       

        const doublereal tauk = tau0 / (1.0 + rho0 * mu[k]);

#pragma omp parallel for
            for (integer i = 0; i < n; ++i) 
            {                
                x[i] += tauk * r[i];
            }
        
    }

    delete[] mu;

    // calculate delta
    residual_for_Cheb(val, col_ind, row_ptr, b, x, r, n);

    doublereal delta = NormaV(r, n) / r_start;

    delete[] r;
    


    // Итерационное уточнение нижней границы спектра.
    // 1.
    doublereal etta = lo / hi;
    doublereal rho1 = (1.0 - sqrt(etta)) / (1.0 + sqrt(etta));
    doublereal qp = 2.0 * pow(rho1, 1.0 * degree) / (1.0 + pow(rho1, 2.0 * degree));
    // 2.
    doublereal y1 = delta / qp;
    doublereal y2 = log(y1 + sqrt(y1 * y1 - 1.0));
    // 3.
    doublereal x_zv = cosh(y2 / (1.0 * degree));
    // 4.
    doublereal etta_new = 0.5 * (1.0 + etta) - 0.5 * (1.0 - etta) * x_zv;
    // 5. 
    //std::cout << "x_zv=" << x_zv << " y2=" << y2 << " y1=" << y1 << " qp=" << qp << " rho1=" << rho1 << " etta=" << etta << " delta=" << delta << std::endl;
    //std::cout << "apostoriory lo = " << lo << "  etta_new=" << etta_new << "  level=" << level << std::endl;
   // getchar();

    if (etta_new > 1.0e-20) {

        switch (iVar) {
        case TEMP:
            level_Chebyshev_info[level].lo_Temp = etta_new * level_Chebyshev_info[level].hi_Temp;
            break;
        case PAM:
            level_Chebyshev_info[level].lo_Pressure = etta_new * level_Chebyshev_info[level].hi_Pressure;
            break;
        case TOTALDEFORMATIONVAR:
            level_Chebyshev_info[level].lo_Stress = etta_new * level_Chebyshev_info[level].hi_Stress;
            break;
        case VELOCITY_X_COMPONENT:
            level_Chebyshev_info[level].lo_Vx = etta_new * level_Chebyshev_info[level].hi_Vx;
            break;
        case VELOCITY_Y_COMPONENT:
            level_Chebyshev_info[level].lo_Vy = etta_new * level_Chebyshev_info[level].hi_Vy;
            break;
        case VELOCITY_Z_COMPONENT:
            level_Chebyshev_info[level].lo_Vz = etta_new * level_Chebyshev_info[level].hi_Vz;
            break;
        case NUSHA:
            level_Chebyshev_info[level].lo_nu = etta_new * level_Chebyshev_info[level].hi_nu;
            break;
        case TURBULENT_KINETIK_ENERGY:
            level_Chebyshev_info[level].lo_k = etta_new * level_Chebyshev_info[level].hi_k;
            break;
        case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
            level_Chebyshev_info[level].lo_omega = etta_new * level_Chebyshev_info[level].hi_omega;
            break;
        case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
            level_Chebyshev_info[level].lo_k_for_ke = etta_new * level_Chebyshev_info[level].hi_k_for_ke;
            break;
        case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
            level_Chebyshev_info[level].lo_epsilon = etta_new * level_Chebyshev_info[level].hi_epsilon;
            break;
        case RE_THETA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_ReTheta = etta_new * level_Chebyshev_info[level].hi_ReTheta;
            break;
        case GAMMA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_gamma = etta_new * level_Chebyshev_info[level].hi_gamma;
            break;
        default:
            level_Chebyshev_info[level].lo_Speed = etta_new * level_Chebyshev_info[level].hi_Speed;
            break;
        }
    }

}

// Вычисление невязки (специально для метода Чебышева).
inline void residual_for_Cheb(Ak1*& Amat, integer istartq, integer iendq,  doublereal*& b, doublereal*& x, doublereal*& r, 
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{


    const integer startpos = istartq + iadd;
    const integer endpos = iendq + iadd;

    //MatrixCRSByVector(val, col_ind, row_ptr, x, r, n);

#pragma omp parallel for 
    for (integer i = startpos; i <= endpos; ++i) {

       // const integer istr = i - iadd;

        const integer is1 = row_ptr_start[i];// +1;
        const integer is2 = row_ptr_end[i];

        doublereal sum = 0.0;

        for (integer j = is1; j <= is2; ++j)
        {
            if (Amat[j].j != Amat[j].i) {
                sum += Amat[j].aij * x[Amat[j].j];
            }
            else {
                sum += (1.0/Amat[j].aij) * x[Amat[j].j];
            }
        }
        r[i-startpos] = sum;
    }

#pragma omp parallel for
    for (integer i = startpos; i <= endpos; ++i)
    {
        const integer istr = i - iadd;
        const integer imover = i - startpos;
        r[imover] = b[istr] - r[imover];

    }

} // residual_for_Cheb

// нумерация  x и ap начинается с нуля.
// Вычисление невязки (специально для метода Чебышева).
inline void Matrix_by_vector_for_Cheb(Ak2& Amat, integer istartq, integer iendq, doublereal*& x, doublereal*& ap,
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer level)
{


    const integer startpos = istartq + iadd;
    const integer endpos = iendq + iadd;

    //MatrixCRSByVector(val, col_ind, row_ptr, x, r, n);

#pragma omp parallel for 
    for (integer i = startpos; i <= endpos; ++i) {

        // const integer istr = i - iadd;

        const integer is1 = row_ptr_start[i];// +1;
        const integer is2 = row_ptr_end[i];

        doublereal sum = (1.0 / Amat.aij[is1]) * x[Amat.j[is1]-1];


        for (integer j = is1 + 1; j <= is2; ++j)
        {
            sum += Amat.aij[j] * x[Amat.j[j]-1];
        }
        ap[i-1] = sum;
    }

} // Matrix_by_vector_for_Cheb

// Вычисление невязки (специально для метода Чебышева).
inline void residual_for_Cheb(Ak2 &Amat, integer istartq, integer iendq, doublereal*& b, doublereal*& x, doublereal*& r,
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer level)
{


    const integer startpos = istartq + iadd;
    const integer endpos = iendq + iadd;

    //MatrixCRSByVector(val, col_ind, row_ptr, x, r, n);

#pragma omp parallel for 
    for (integer i = startpos; i <= endpos; ++i) {

       // const integer istr = i - iadd;

        const integer is1 = row_ptr_start[i];// +1;
        const integer is2 = row_ptr_end[i];

        doublereal sum = (1.0 / Amat.aij[is1]) * x[Amat.j[is1]];
       

        for (integer j = is1+1; j <= is2; ++j)
        {            
             sum += Amat.aij[j] * x[Amat.j[j]];
        }
        r[i - startpos] = b[i - iadd] - sum;
    }

} // residual_for_Cheb

// Вычисление невязки (специально для метода Чебышева).
inline doublereal normaV_Cheb(Ak2& Amat, integer istartq, integer iendq, doublereal*& b, doublereal*& x, 
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer level)
{
    const integer startpos = istartq + iadd;
    const integer endpos = iendq + iadd;

    //MatrixCRSByVector(val, col_ind, row_ptr, x, r, n);

    doublereal sum1 = 0.0;

#pragma omp parallel for reduction(+: sum1)
    for (integer i = startpos; i <= endpos; ++i) {

        // const integer istr = i - iadd;

        const integer is1 = row_ptr_start[i];// +1;
        const integer is2 = row_ptr_end[i];

        doublereal sum = (1.0 / Amat.aij[is1]) * x[Amat.j[is1]];
        

        for (integer j = is1 + 1; j <= is2; ++j)
        {
            sum += Amat.aij[j] * x[Amat.j[j]];
        }
        doublereal v = b[i - iadd] - sum;
        sum1 += v * v;
    }

    return (sqrt(sum1/(1.0*(iendq - istartq + 1))));

} // normaV_Cheb

// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.
// Ak1 не тестировался!!!
doublereal spectral_radius(Ak1*& Amat, integer istartq, integer iendq, 
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool scale) {
    //std::cout << "spectral radius" << std::endl;


    doublereal radius;

    //if (power_iters <= 0)
    {
        // Use Gershgorin disk theorem.
        radius = 0;

#pragma omp parallel
        {
            doublereal emax = 0;
            doublereal  dia = 1.0;

            const integer startpos = istartq + iadd;
            const integer endpos = iendq + iadd;

#pragma omp for nowait
            for (integer i = startpos; i <= endpos; ++i) {
                doublereal s = 0;

                for (integer j = row_ptr_start[i], e = row_ptr_end[i]; j <= e; ++j) {
                    integer  c = Amat[j].j;
                    doublereal v = Amat[j].aij;

                    if (c == Amat[j].i) v = 1 / v;

                    s += fabs(v);// норма числа

                    if (scale && c == i) dia = v;
                }

                if (scale) s *= 1.0 / (fabs(dia));

                emax = std::max(emax, s);
            }

#pragma omp critical
            radius = std::max(radius, emax);
        }
    }
    //std::cout << "spectral radius" << std::endl;

    return radius < 0 ? static_cast<doublereal>(2.0) : radius;
}


// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.

doublereal spectral_radius(Ak2 &Amat, integer istartq, integer iendq,
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool scale, integer level) {
    //std::cout << "spectral radius" << std::endl;


    doublereal radius;

    //if (power_iters <= 0)
    {
        // Use Gershgorin disk theorem.
        radius = 0;

//#pragma omp parallel
        {
            doublereal emax = 0;
            doublereal  dia = 1.0;

            const integer startpos = istartq + iadd;
            const integer endpos = iendq + iadd;

//#pragma omp for nowait
            for (integer i = startpos; i <= endpos; ++i) {
                doublereal s = 0;


                 

                for (integer j = row_ptr_start[i], e = row_ptr_end[i]; j <= e; ++j) {
                    //integer  c = Amat.j[j];
                    doublereal v = Amat.aij[j];

                    //if (level==1) {
                       //std::cout << "ilev="<< level << "  jstart = " << row_ptr_start[i] << " i = " << i << " " << Amat.aij[j] << "j = " << Amat.j[j] << std::endl;
                    //}

                    //if (c == Amat.i[j]) v = 1 / v;
                    if (j == row_ptr_start[i]) v = 1 / v;// первый элемент в строке диагональ.

                    s += fabs(v);// норма числа

                    if (scale && j == row_ptr_start[i]) dia = v;
                }
               // if (level == 1) {
                    //getchar();
                //}

                if (scale) s *= 1.0 / (fabs(dia));

                emax = std::max(emax, s);
            }

//#pragma omp critical
            radius = std::max(radius, emax);
        }
    }
    //std::cout << "spectral radius" << std::endl;

    return radius < 0 ? static_cast<doublereal>(2.0) : radius;
}

// 26.01.2022
template <typename doublerealT>
void chebyshev(Ak1*& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b,
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer iVar, integer level)
{
    // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.


   /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5; 64;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
    bool scale;


    higher = 1.0f;
    lower = 1.0f / 30;
   // power_iters = 0;
    scale = false;


    doublereal hi, lo;



    // Спектральный радиус вычисляем если он еще не был вычислен ранее по теореме Гершгорина.
    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Temp) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Temp = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Temp = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Temp = true;
        }
        hi = level_Chebyshev_info[level].hi_Temp;
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Pressure) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Pressure = true;
        }
        hi = level_Chebyshev_info[level].hi_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Stress) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Stress = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Stress = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Stress = true;
        }
        hi = level_Chebyshev_info[level].hi_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vx) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vx = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Vx = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vx = true;
        }
        hi = level_Chebyshev_info[level].hi_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vy) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vy = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Vy = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vy = true;
        }
        hi = level_Chebyshev_info[level].hi_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vz) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vz = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Vz = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vz = true;
        }
        hi = level_Chebyshev_info[level].hi_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_nu) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_nu = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_nu = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_nu = true;
        }
        hi = level_Chebyshev_info[level].hi_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_k = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k = true;
        }
        hi = level_Chebyshev_info[level].hi_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_omega) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_omega = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_omega = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_omega = true;
        }
        hi = level_Chebyshev_info[level].hi_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k_for_ke) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k_for_ke = true;
        }
        hi = level_Chebyshev_info[level].hi_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_epsilon) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_epsilon = true;
        }
        hi = level_Chebyshev_info[level].hi_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_ReTheta) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_ReTheta = true;
        }
        hi = level_Chebyshev_info[level].hi_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_gamma) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_gamma = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_gamma = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_gamma = true;
        }
        hi = level_Chebyshev_info[level].hi_gamma;
        break;
    default:

        if (!level_Chebyshev_info[level].bGershgorin_calc_Speed) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Speed = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true);
            }
            else {
                level_Chebyshev_info[level].hi_Speed = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Speed = true;
        }
        hi = level_Chebyshev_info[level].hi_Speed;
        break;
    }

    integer n = iendq - istartq + 1;

    doublereal* r = new doublereal[n];



    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bstart_lo_Temp) {
            level_Chebyshev_info[level].lo_Temp = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Temp = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Temp = true;
        lo = level_Chebyshev_info[level].lo_Temp;
        //if (lo < 1.0e-14)
        {
            //std::cout << "incomming lo = " << lo << " level=" << level << std::endl;
          //  getchar();
        }
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bstart_lo_Pressure) {
            level_Chebyshev_info[level].lo_Pressure = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Pressure = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Pressure = true;
        lo = level_Chebyshev_info[level].lo_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bstart_lo_Stress) {
            level_Chebyshev_info[level].lo_Stress = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Stress = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Stress = true;
        lo = level_Chebyshev_info[level].lo_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vx) {
            level_Chebyshev_info[level].lo_Vx = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vx = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vx = true;
        lo = level_Chebyshev_info[level].lo_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vy) {
            level_Chebyshev_info[level].lo_Vy = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vy = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vy = true;
        lo = level_Chebyshev_info[level].lo_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vz) {
            level_Chebyshev_info[level].lo_Vz = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vz = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vz = true;
        lo = level_Chebyshev_info[level].lo_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bstart_lo_nu) {
            level_Chebyshev_info[level].lo_nu = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_nu = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_nu = true;
        lo = level_Chebyshev_info[level].lo_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bstart_lo_k) {
            level_Chebyshev_info[level].lo_k = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k = true;
        lo = level_Chebyshev_info[level].lo_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bstart_lo_omega) {
            level_Chebyshev_info[level].lo_omega = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_omega = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_omega = true;
        lo = level_Chebyshev_info[level].lo_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_k_for_ke) {
            level_Chebyshev_info[level].lo_k_for_ke = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k_for_ke = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k_for_ke = true;
        lo = level_Chebyshev_info[level].lo_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_epsilon) {
            level_Chebyshev_info[level].lo_epsilon = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_epsilon = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_epsilon = true;
        lo = level_Chebyshev_info[level].lo_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_ReTheta) {
            level_Chebyshev_info[level].lo_ReTheta = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_ReTheta = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_ReTheta = true;
        lo = level_Chebyshev_info[level].lo_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_gamma) {
            level_Chebyshev_info[level].lo_gamma = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_gamma = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_gamma = true;
        lo = level_Chebyshev_info[level].lo_gamma;
        break;
    default:
        if (!level_Chebyshev_info[level].bstart_lo_Speed) {
            level_Chebyshev_info[level].lo_Speed = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Speed = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Speed = true;
        lo = level_Chebyshev_info[level].lo_Speed;
        break;
    }



    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();



    doublereal r_start;

    const doublereal tau0 = 2.0 / (hi + lo);
    const doublereal ksi = lo / hi;
    const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

   
    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);


    const integer startpos = istartq + iadd;
    const integer endpos = iendq + iadd;

    for (unsigned k = 0; k < degree; ++k) {
        //backend::residual(b, A, x, *r);
        residual_for_Cheb(Amat, istartq, iendq, b, x, r, row_ptr_start, row_ptr_end, iadd);

        if (k == 0) r_start = NormaV(r, n);

        //if (prm.scale) backend::vmul(one, *M, *r, zero, *r);       

        const doublereal tauk = tau0 / (1.0 + rho0 * mu[k]);

//#pragma omp parallel for
  //      for (integer i = 0; i < n; ++i)
    //    {
      //      x[i] += tauk * r[i];
      //  }

#pragma omp parallel for
        for (integer i = startpos; i <= endpos; ++i)
        {
            const integer istr = i - iadd;
            const integer imover = i - startpos;
            
            x[istr] += tauk * r[imover];

        }

    }

    // calculate delta
    residual_for_Cheb(Amat, istartq, iendq, b, x, r, row_ptr_start, row_ptr_end, iadd);

    doublereal delta = NormaV(r, n) / r_start;

    delete[] r;
    delete[] mu;



    // Итерационное уточнение нижней границы спектра.
    // 1.
    doublereal etta = lo / hi;
    doublereal rho1 = (1.0 - sqrt(etta)) / (1.0 + sqrt(etta));
    doublereal qp = 2.0 * pow(rho1, 1.0 * degree) / (1.0 + pow(rho1, 2.0 * degree));
    // 2.
    doublereal y1 = delta / qp;
    doublereal y2 = log(y1 + sqrt(y1 * y1 - 1.0));
    // 3.
    doublereal x_zv = cosh(y2 / (1.0 * degree));
    // 4.
    doublereal etta_new = 0.5 * (1.0 + etta) - 0.5 * (1.0 - etta) * x_zv;
    // 5. 
    //std::cout << "x_zv=" << x_zv << " y2=" << y2 << " y1=" << y1 << " qp=" << qp << " rho1=" << rho1 << " etta=" << etta << " delta=" << delta << std::endl;
    //std::cout << "apostoriory lo = " << lo << "  etta_new=" << etta_new << "  level=" << level << std::endl;
   // getchar();

    if (etta_new > 1.0e-20) {

        switch (iVar) {
        case TEMP:
            level_Chebyshev_info[level].lo_Temp = etta_new * level_Chebyshev_info[level].hi_Temp;
            break;
        case PAM:
            level_Chebyshev_info[level].lo_Pressure = etta_new * level_Chebyshev_info[level].hi_Pressure;
            break;
        case TOTALDEFORMATIONVAR:
            level_Chebyshev_info[level].lo_Stress = etta_new * level_Chebyshev_info[level].hi_Stress;
            break;
        case VELOCITY_X_COMPONENT:
            level_Chebyshev_info[level].lo_Vx = etta_new * level_Chebyshev_info[level].hi_Vx;
            break;
        case VELOCITY_Y_COMPONENT:
            level_Chebyshev_info[level].lo_Vy = etta_new * level_Chebyshev_info[level].hi_Vy;
            break;
        case VELOCITY_Z_COMPONENT:
            level_Chebyshev_info[level].lo_Vz = etta_new * level_Chebyshev_info[level].hi_Vz;
            break;
        case NUSHA:
            level_Chebyshev_info[level].lo_nu = etta_new * level_Chebyshev_info[level].hi_nu;
            break;
        case TURBULENT_KINETIK_ENERGY:
            level_Chebyshev_info[level].lo_k = etta_new * level_Chebyshev_info[level].hi_k;
            break;
        case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
            level_Chebyshev_info[level].lo_omega = etta_new * level_Chebyshev_info[level].hi_omega;
            break;
        case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
            level_Chebyshev_info[level].lo_k_for_ke = etta_new * level_Chebyshev_info[level].hi_k_for_ke;
            break;
        case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
            level_Chebyshev_info[level].lo_epsilon = etta_new * level_Chebyshev_info[level].hi_epsilon;
            break;
        case RE_THETA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_ReTheta = etta_new * level_Chebyshev_info[level].hi_ReTheta;
            break;
        case GAMMA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_gamma = etta_new * level_Chebyshev_info[level].hi_gamma;
            break;
        default:
            level_Chebyshev_info[level].lo_Speed = etta_new * level_Chebyshev_info[level].hi_Speed;
            break;
        }
    }

}


// 26.01.2022
template <typename doublerealT>
void chebyshev(Ak2 &Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b,
    integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer iVar, integer level)
{
    // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.


    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
    bool scale;


    higher = 1.0f;
    lower = 1.0f / 30;
    //lower = 1.0f / 60; // достаточно одной тридцатой.
    //power_iters = 0;
    scale = false;


    doublereal hi, lo;



    // Спектральный радиус вычисляем если он еще не был вычислен ранее по теореме Гершгорина.
    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Temp) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Temp = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Temp = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Temp = true;
        }
        hi = level_Chebyshev_info[level].hi_Temp;
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Pressure) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Pressure = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Pressure = true;
        }
        hi = level_Chebyshev_info[level].hi_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Stress) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Stress = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Stress = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Stress = true;
        }
        hi = level_Chebyshev_info[level].hi_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vx) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vx = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Vx = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vx = true;
        }
        hi = level_Chebyshev_info[level].hi_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vy) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vy = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Vy = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vy = true;
        }
        hi = level_Chebyshev_info[level].hi_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bGershgorin_calc_Vz) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Vz = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Vz = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Vz = true;
        }
        hi = level_Chebyshev_info[level].hi_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_nu) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_nu = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_nu = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_nu = true;
        }
        hi = level_Chebyshev_info[level].hi_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_k = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k = true;
        }
        hi = level_Chebyshev_info[level].hi_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bGershgorin_calc_omega) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_omega = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_omega = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_omega = true;
        }
        hi = level_Chebyshev_info[level].hi_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_k_for_ke) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_k_for_ke = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_k_for_ke = true;
        }
        hi = level_Chebyshev_info[level].hi_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bGershgorin_calc_epsilon) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_epsilon = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_epsilon = true;
        }
        hi = level_Chebyshev_info[level].hi_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_ReTheta) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_ReTheta = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_ReTheta = true;
        }
        hi = level_Chebyshev_info[level].hi_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bGershgorin_calc_gamma) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_gamma = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_gamma = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_gamma = true;
        }
        hi = level_Chebyshev_info[level].hi_gamma;
        break;
    default:

        if (!level_Chebyshev_info[level].bGershgorin_calc_Speed) {
            // Применяем теорему Гершгорина это требует ресурсов.
            if (scale) {
                //M = Backend::copy_vector(diagonal(A, /*invert*/true), backend_prm);
                level_Chebyshev_info[level].hi_Speed = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, true, level);
            }
            else {
                level_Chebyshev_info[level].hi_Speed = spectral_radius(Amat, istartq, iendq, row_ptr_start, row_ptr_end, iadd, false, level);
            }
            level_Chebyshev_info[level].bGershgorin_calc_Speed = true;
        }
        hi = level_Chebyshev_info[level].hi_Speed;
        break;
    }

    integer n = iendq - istartq + 1;

    doublereal* r = new doublereal[n];



    switch (iVar) {
    case TEMP:
        if (!level_Chebyshev_info[level].bstart_lo_Temp) {
            level_Chebyshev_info[level].lo_Temp = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Temp = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Temp = true;
        lo = level_Chebyshev_info[level].lo_Temp;
        //if (lo < 1.0e-14)
        {
            //std::cout << "incomming lo = " << lo << " level=" << level << std::endl;
          //  getchar();
        }
        break;
    case PAM:
        if (!level_Chebyshev_info[level].bstart_lo_Pressure) {
            level_Chebyshev_info[level].lo_Pressure = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Pressure = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Pressure = true;
        lo = level_Chebyshev_info[level].lo_Pressure;
        break;
    case TOTALDEFORMATIONVAR:
        if (!level_Chebyshev_info[level].bstart_lo_Stress) {
            level_Chebyshev_info[level].lo_Stress = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Stress = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Stress = true;
        lo = level_Chebyshev_info[level].lo_Stress;
        break;
    case VELOCITY_X_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vx) {
            level_Chebyshev_info[level].lo_Vx = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vx = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vx = true;
        lo = level_Chebyshev_info[level].lo_Vx;
        break;
    case VELOCITY_Y_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vy) {
            level_Chebyshev_info[level].lo_Vy = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vy = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vy = true;
        lo = level_Chebyshev_info[level].lo_Vy;
        break;
    case VELOCITY_Z_COMPONENT:
        if (!level_Chebyshev_info[level].bstart_lo_Vz) {
            level_Chebyshev_info[level].lo_Vz = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_Vz = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Vz = true;
        lo = level_Chebyshev_info[level].lo_Vz;
        break;
    case NUSHA:
        if (!level_Chebyshev_info[level].bstart_lo_nu) {
            level_Chebyshev_info[level].lo_nu = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_nu = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_nu = true;
        lo = level_Chebyshev_info[level].lo_nu;
        break;
    case TURBULENT_KINETIK_ENERGY:
        if (!level_Chebyshev_info[level].bstart_lo_k) {
            level_Chebyshev_info[level].lo_k = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k = true;
        lo = level_Chebyshev_info[level].lo_k;
        break;
    case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
        if (!level_Chebyshev_info[level].bstart_lo_omega) {
            level_Chebyshev_info[level].lo_omega = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_omega = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_omega = true;
        lo = level_Chebyshev_info[level].lo_omega;
        break;
    case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_k_for_ke) {
            level_Chebyshev_info[level].lo_k_for_ke = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_k_for_ke = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_k_for_ke = true;
        lo = level_Chebyshev_info[level].lo_k_for_ke;
        break;
    case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
        if (!level_Chebyshev_info[level].bstart_lo_epsilon) {
            level_Chebyshev_info[level].lo_epsilon = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_epsilon = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_epsilon = true;
        lo = level_Chebyshev_info[level].lo_epsilon;
        break;
    case RE_THETA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_ReTheta) {
            level_Chebyshev_info[level].lo_ReTheta = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_ReTheta = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_ReTheta = true;
        lo = level_Chebyshev_info[level].lo_ReTheta;
        break;
    case GAMMA_LANGTRY_MENTER:
        if (!level_Chebyshev_info[level].bstart_lo_gamma) {
            level_Chebyshev_info[level].lo_gamma = hi * lower;

            // Соотношение Релея - Ритца.
                    //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
                    //level_Chebyshev_info[level].lo_gamma = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_gamma = true;
        lo = level_Chebyshev_info[level].lo_gamma;
        break;
    default:
        if (!level_Chebyshev_info[level].bstart_lo_Speed) {
            level_Chebyshev_info[level].lo_Speed = hi * lower;

            // Соотношение Релея - Ритца.
            //MatrixCRSByVector(val, col_ind, row_ptr, b, r, n);
            //level_Chebyshev_info[level].lo_Speed = Scal(r, b, n) / Scal(b, b, n);
        }
        level_Chebyshev_info[level].bstart_lo_Speed = true;
        lo = level_Chebyshev_info[level].lo_Speed;
        break;
    }



    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();



    doublereal r_start;

    const doublereal tau0 = 2.0 / (hi + lo);
    const doublereal ksi = lo / hi;
    const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

   

     // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);
    


    const integer startpos = istartq + iadd;
    const integer endpos = iendq + iadd;

    for (unsigned k = 0; k < degree; ++k) {
        //backend::residual(b, A, x, *r);
        residual_for_Cheb(Amat, istartq, iendq, b, x, r, row_ptr_start, row_ptr_end, iadd, level);

        if (k == 0) r_start = NormaV(r, n);

        //if (prm.scale) backend::vmul(one, *M, *r, zero, *r);       

        const doublereal tauk = tau0 / (1.0 + rho0 * mu[k]);

        //#pragma omp parallel for
          //      for (integer i = 0; i < n; ++i)
            //    {
              //      x[i] += tauk * r[i];
              //  }

#pragma omp parallel for
        for (integer i = startpos; i <= endpos; ++i)
        {
            const integer istr = i - iadd;
            const integer imover = i - startpos;

            x[istr] += tauk * r[imover];

        }

    }

    delete[] r;
    delete[] mu;


    // calculate delta
    doublereal delta = normaV_Cheb(Amat, istartq, iendq, b, x, row_ptr_start, row_ptr_end, iadd, level) / r_start;

   



    // Итерационное уточнение нижней границы спектра.
    // 1.
    doublereal etta = lo / hi;
    doublereal rho1 = (1.0 - sqrt(etta)) / (1.0 + sqrt(etta));
    doublereal qp = 2.0 * pow(rho1, 1.0 * degree) / (1.0 + pow(rho1, 2.0 * degree));
    // 2.
    doublereal y1 = delta / qp;
    doublereal y2 = log(y1 + sqrt(y1 * y1 - 1.0));
    // 3.
    doublereal x_zv = cosh(y2 / (1.0 * degree));
    // 4.
    doublereal etta_new = 0.5 * (1.0 + etta) - 0.5 * (1.0 - etta) * x_zv;
    // 5. 
    //std::cout << "x_zv=" << x_zv << " y2=" << y2 << " y1=" << y1 << " qp=" << qp << " rho1=" << rho1 << " etta=" << etta << " delta=" << delta << std::endl;
    //std::cout << "apostoriory lo = " << lo << "  etta_new=" << etta_new << "  level=" << level << std::endl;
   // getchar();

    if (etta_new > 1.0e-20) {

        switch (iVar) {
        case TEMP:
            level_Chebyshev_info[level].lo_Temp = etta_new * level_Chebyshev_info[level].hi_Temp;
            break;
        case PAM:
            level_Chebyshev_info[level].lo_Pressure = etta_new * level_Chebyshev_info[level].hi_Pressure;
            break;
        case TOTALDEFORMATIONVAR:
            level_Chebyshev_info[level].lo_Stress = etta_new * level_Chebyshev_info[level].hi_Stress;
            break;
        case VELOCITY_X_COMPONENT:
            level_Chebyshev_info[level].lo_Vx = etta_new * level_Chebyshev_info[level].hi_Vx;
            break;
        case VELOCITY_Y_COMPONENT:
            level_Chebyshev_info[level].lo_Vy = etta_new * level_Chebyshev_info[level].hi_Vy;
            break;
        case VELOCITY_Z_COMPONENT:
            level_Chebyshev_info[level].lo_Vz = etta_new * level_Chebyshev_info[level].hi_Vz;
            break;
        case NUSHA:
            level_Chebyshev_info[level].lo_nu = etta_new * level_Chebyshev_info[level].hi_nu;
            break;
        case TURBULENT_KINETIK_ENERGY:
            level_Chebyshev_info[level].lo_k = etta_new * level_Chebyshev_info[level].hi_k;
            break;
        case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
            level_Chebyshev_info[level].lo_omega = etta_new * level_Chebyshev_info[level].hi_omega;
            break;
        case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
            level_Chebyshev_info[level].lo_k_for_ke = etta_new * level_Chebyshev_info[level].hi_k_for_ke;
            break;
        case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
            level_Chebyshev_info[level].lo_epsilon = etta_new * level_Chebyshev_info[level].hi_epsilon;
            break;
        case RE_THETA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_ReTheta = etta_new * level_Chebyshev_info[level].hi_ReTheta;
            break;
        case GAMMA_LANGTRY_MENTER:
            level_Chebyshev_info[level].lo_gamma = etta_new * level_Chebyshev_info[level].hi_gamma;
            break;
        default:
            level_Chebyshev_info[level].lo_Speed = etta_new * level_Chebyshev_info[level].hi_Speed;
            break;
        }
    }

}




// применяется для уравнения поправки давления. 
// Работает только на структурированной сетке данная версия.
// 05.02.2022.
// 
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
void chebyshev3Dnow(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
    doublereal*& x, integer maxelm, integer maxbound, integer iVar,
    bool bprintMessage, bool bsolidT) {
    // К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
    // скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
    //printf("chebyshev3Dnow incomming...\n"); // debug.
    //system("pause");

    // Параметр релаксации лежит в интервале от 0.0 до 2.0.
    // Верхняя релаксация способна существенно ускорить вычислительный процесс. 
    // Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
    // rURF=1.5 можно использовать в уравнении на поправку давления.
    // В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
    // к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
    // в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
    // можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
    // коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
    // поэтому наверное лучше не вводить коээффициент 
    // верхней релаксации равный 1.5 в уравнениях на скорость. 

     // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

#ifdef _OPENMP 
    omp_set_num_threads(number_cores()); // установка числа потоков
#endif

    if (bsolidT) {
        // диагональное предобуславливание
        // для температуры в твёрдом теле.
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae /= sl[i].ap;
            sl[i].aw /= sl[i].ap;
            sl[i].an /= sl[i].ap;
            sl[i].as /= sl[i].ap;
            sl[i].at /= sl[i].ap;
            sl[i].ab /= sl[i].ap;

            b[sl[i].iP] /= sl[i].ap;
            sl[i].ap = 1.0;
        }

        for (int i = 0; i < maxbound; ++i) {
            slb[i].ai /= slb[i].aw;
            b[slb[i].iW] /= slb[i].aw;

            slb[i].aw = 1.0;
        }
    }

    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
   // bool scale;


    higher = 1.0f;
    lower = 1.0f / 3.0f;// 1.0f / 30;
    //if (inumiterSIMPLE371 > 30) {
      //  lower = free_debug_parametr1;
   // }
    //lower = 1.0f / 60; // достаточно одной тридцатой.
    //power_iters = 0;
   // scale = false;


    doublereal hi, lo;

    hi = spectral_radius(sl, slb,  maxelm,  maxbound);

    lo = hi * lower;


    doublereal rURF = 1.0; // параметр верхней релаксации
    
    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();



    



    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);

    // пороговое значение невязки
    doublereal eps = 1e-40;
    
     
    integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.


    const integer size0 = maxelm + maxbound;

    doublereal* x_new = new doublereal[size0];
#pragma omp parallel for
    for (integer i = 0; i < size0; ++i)
    {
        x_new[i]=0.0;
    }

    {
        // Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
        // опираемся на опыт теплопередачи.
        doublereal c1 = 0.026; // воздух
        doublereal c2 = 2000; // алмаз

        // см. Николаев
        if (bprintMessage) {
            std::cout << " number equations " << maxelm + maxbound << std::endl;
        }
        kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
        if (bprintMessage) {
            std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
        }
    }

  

    // Улучшает попадание в кеш ускорение на 12.5%.
    doublereal* data = new doublereal[maxelm*7];
    int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
    for (int i = 0; i < maxelm; ++i) {

        const int is = i * 7;
        data[is] = sl[i].ae;
        data[is + 1] = sl[i].aw;
        data[is + 2] = sl[i].an;
        data[is + 3] = sl[i].as;
        data[is + 4] = sl[i].at;
        data[is + 5] = sl[i].ab;
        data[is + 6] = sl[i].ap;

        col_ind[is] = sl[i].iE;
        col_ind[is + 1] = sl[i].iW;
        col_ind[is + 2] = sl[i].iN;
        col_ind[is + 3] = sl[i].iS;
        col_ind[is + 4] = sl[i].iT;
        col_ind[is + 5] = sl[i].iB;
        col_ind[is + 6] = sl[i].iP;
    }


    doublereal* data_b = new doublereal[maxbound * 2];
    int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
    for (int i = 0; i < maxbound; ++i) {

        const int is = i * 2;
        data_b[is] = slb[i].ai;
        data_b[is + 1] = slb[i].aw;

        if (slb[i].iI > -1) {
            col_ind_b[is] = slb[i].iI;
        }
        else {
            data_b[is] = 0.0;
            col_ind_b[is] = slb[i].iW;
        }
        col_ind_b[is + 1] = slb[i].iW;
    }

    doublereal dmax = 1.0, dmax2=0.0, delta=0.0;
    integer iout = 0;
    while ((dmax > eps) && (j < kend*degree)) {
       
        if (bsolidT) {
            if (dmax < 1.0e-4) break; // Вычисление сошлось.
        }

        bool bmem = true;

        doublereal r_start;

        const doublereal tau0 = 2.0 / (hi + lo);
        const doublereal ksi = lo / hi;
        const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

        for (int k = 0; k < degree; ++k) {

           

            rURF= tau0 / (1.0 + rho0 * mu[k]);

            dmax = 0.0;
            dmax2 = 0.0;
            delta = 0.0;

            doublereal dmax_loc = 0.0, dmax2_loc=0.0; 
            
            if ((k == 0) || (k == degree - 1)) {

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
                for (int i = 0; i < maxelm; ++i) {



                    //long double residualQ = (sl[i].ae * x[sl[i].iE] +
                      //              sl[i].aw * x[sl[i].iW] +
                        //            sl[i].an * x[sl[i].iN] +
                          //          sl[i].as * x[sl[i].iS] +
                            //        sl[i].at * x[sl[i].iT] +
                              //      sl[i].ab * x[sl[i].iB] +
                                //    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

                    const int is = i * 7;
                    const int is1 = is + 1;
                    const int is2 = is + 2;
                    const int is3 = is + 3;
                    const int is4 = is + 4;
                    const int is5 = is + 5;
                    const int is6 = is + 6;
                    const int it6 = col_ind[is6];

                    const long double residualQ = (data[is] * x[col_ind[is]] +
                        data[is1] * x[col_ind[is1]] +
                        data[is2] * x[col_ind[is2]] +
                        data[is3] * x[col_ind[is3]] +
                        data[is4] * x[col_ind[is4]] +
                        data[is5] * x[col_ind[is5]] +
                        b[it6]) - data[is6] * x[it6];

                    const doublereal residual = static_cast<doublereal>(residualQ);

                    // Чебышева
                    dmax_loc += fabs((residual)); // сумма модулей
                    dmax2_loc += residual * residual;
                    x_new[it6] = x[it6] + rURF * (residual);


                }
            }
            else {

                // Редукция медленная и здесь мы от неё отказываемся в пользу скорости
                // вычисления т.к. внутри она ненужна.

#pragma omp parallel for 
                for (int i = 0; i < maxelm; ++i) {



                    //long double residualQ = (sl[i].ae * x[sl[i].iE] +
                      //              sl[i].aw * x[sl[i].iW] +
                        //            sl[i].an * x[sl[i].iN] +
                          //          sl[i].as * x[sl[i].iS] +
                            //        sl[i].at * x[sl[i].iT] +
                              //      sl[i].ab * x[sl[i].iB] +
                                //    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

                    const int is = i * 7;
                    const int is1 = is + 1;
                    const int is2 = is + 2;
                    const int is3 = is + 3;
                    const int is4 = is + 4;
                    const int is5 = is + 5;
                    const int is6 = is + 6;
                    const int it6 = col_ind[is6];

                    const long double residualQ = (data[is] * x[col_ind[is]] +
                        data[is1] * x[col_ind[is1]] +
                        data[is2] * x[col_ind[is2]] +
                        data[is3] * x[col_ind[is3]] +
                        data[is4] * x[col_ind[is4]] +
                        data[is5] * x[col_ind[is5]] +
                        b[it6]) - data[is6] * x[it6];

                    const doublereal residual = static_cast<doublereal>(residualQ);

                    
                    x_new[it6] = x[it6] + rURF * (residual);

                }

            }



            dmax += dmax_loc;
            dmax2 += dmax2_loc;

            dmax_loc = 0.0;
            dmax2_loc = 0.0;


            if ((k == 0) || (k == degree - 1)) {

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
                for (int i = 0; i < maxbound; ++i) {


                    const int is = i * 2;
                    const int it1 = col_ind_b[is];//iI
                    const int it2 = col_ind_b[is + 1];//iW                   


                    const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
                    // Чебышева
                    dmax_loc += fabs(residual); // сумма модулей.
                    dmax2_loc += residual * residual;

                    if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
                    else x_new[it2] = x[it2] + rURF * (residual);


                }

            }
            else {
                // Редукция медленная и здесь мы от неё отказываемся в пользу скорости
                // вычисления т.к. внутри она ненужна.

#pragma omp parallel for
                for (int i = 0; i < maxbound; ++i) {


                    const int is = i * 2;
                    const int it1 = col_ind_b[is];//iI
                    const int it2 = col_ind_b[is + 1];//iW                   


                    const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];

                    if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
                    else x_new[it2] = x[it2] + rURF * (residual);


                }
            }



            dmax += dmax_loc;
            dmax2 += dmax2_loc;

            //*std::min(1.0, 17000.0 / size0) 
            if (bsolidT) {
                if (iout == 0) {
                    eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
                    if (j == 0) {
                        std::cout << "start residual =" << dmax << std::endl;
                    }
                }
            }
            else {
                if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
            }

            if (bprintMessage) {
                if ((iVar == PAM)||(iVar == TEMP)) {
                    if ((k == 0) || (k == degree - 1)) {
                        //dmax/=maxelm;
                        if (j % degree == 0) {
#if doubleintprecision == 1
                            printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
                            printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

                        }
                    }
                }
            }

#pragma omp parallel for
            for (integer i = 0; i < size0; ++i)
            {
                x[i] = x_new[i];
            }

            if ((bmem)&&(k == 0)) r_start = sqrt(dmax2 / (size0));
            if (k == degree-1) delta= sqrt(dmax2 / (size0));

            if (bprintMessage) {
                if ((k == 0) || (k == degree - 1)) {
                    if (j % degree == 0) {
                        std::cout << std::endl;
                    }
                }
            }

            j++;
        }// degree chebyshev


        delta /= r_start;

       
        // Корректировка нижней границы спектра lo линейного оператора А.
        // Жуков, Феодоритова. 13.02.2022
        correct_lo(lo, bmem, degree, hi, delta, bprintMessage);
        

        //getchar();
        
        ++iout;
    }


    delete[] data;
    delete[] col_ind;

    delete[] data_b;
    delete[] col_ind_b;

    std::cout << " " << lo << " " << j << " ";

    if (j > kend * degree - 100) {
        std::cout << "\n chebyshev divergence " << dmax << " ";
        system("pause");
    }

    delete[] x_new;

    if (bprintMessage) {
        if ((iVar == PAM) || (iVar == TEMP)) {
            printf("calc complete...\n");
            system("pause");
        }
    }
    //printf("4000 %e \n", dmax);
    //system("pause");

} // chebyshev3Dnow








// Метод тяжёлого шарика не работает. По видимому нарушает закон сохранения 11.02.2022.. 
// Сделан на структурированной сетке. 
// 05.02.2022.
// 
// Метод сопряжённых градиентов предобусловленный методом Чебышева. 12.02.2022.
void cg_chebyshev_3Dnow(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
    doublereal*& x, integer maxelm, integer maxbound, integer iVar,
    bool bprintMessage, bool bsolidT) {
    // К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
    // скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
    //printf("chebyshev3Dnow incomming...\n"); // debug.
    //system("pause");

    // Параметр релаксации лежит в интервале от 0.0 до 2.0.
    // Верхняя релаксация способна существенно ускорить вычислительный процесс. 
    // Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
    // rURF=1.5 можно использовать в уравнении на поправку давления.
    // В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
    // к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
    // в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
    // можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
    // коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
    // поэтому наверное лучше не вводить коээффициент 
    // верхней релаксации равный 1.5 в уравнениях на скорость. 

     // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

#ifdef _OPENMP 
    omp_set_num_threads(number_cores()); // установка числа потоков
#endif

    
    // Чистый cg лучше сходится без деления на диагональ.
    if (0&&bsolidT) {
        // диагональное предобуславливание
        // для температуры в твёрдом теле.
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae /= sl[i].ap;
            sl[i].aw /= sl[i].ap;
            sl[i].an /= sl[i].ap;
            sl[i].as /= sl[i].ap;
            sl[i].at /= sl[i].ap;
            sl[i].ab /= sl[i].ap;

            b[sl[i].iP] /= sl[i].ap;
            sl[i].ap = 1.0;
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            slb[i].ai /= slb[i].aw;
            if (slb[i].iW < maxelm) {

                // Граничный узел индексируется как внутренний узел !!! фатальная ошибка.
                std::cout << "fatal error iW bound < maxelm. In function cg_chebyshev_3Dnow.\n";
                system("PAUSE");

            }
            b[slb[i].iW] /= slb[i].aw;

            slb[i].aw = 1.0;
        }
    }

    std::vector<double> D(maxelm + maxbound, 1.0);

    if (bsolidT) {

        // Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            D[sl[i].iP] = 1 / sqrt(fabs(sl[i].ap));
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            D[slb[i].iW] = 1 / sqrt(fabs(slb[i].aw));
        }

        // Then, apply the scaling in-place:
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae *= D[sl[i].iP] * D[sl[i].iE];
            sl[i].aw *= D[sl[i].iP] * D[sl[i].iW];
            sl[i].an *= D[sl[i].iP] * D[sl[i].iN];
            sl[i].as *= D[sl[i].iP] * D[sl[i].iS];
            sl[i].at *= D[sl[i].iP] * D[sl[i].iT];
            sl[i].ab *= D[sl[i].iP] * D[sl[i].iB];
            sl[i].ap *= D[sl[i].iP] * D[sl[i].iP];

            b[sl[i].iP] *= D[sl[i].iP];
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            if (slb[i].iI > -1) {
                slb[i].ai *= D[slb[i].iW] * D[slb[i].iI];
            }
            slb[i].aw *= D[slb[i].iW] * D[slb[i].iW];

            b[slb[i].iW] *= D[slb[i].iW];
        }

        
       
    }
    
    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
   // bool scale;


    higher = 1.0f;
    lower = 1.0f / 3.0f;// 1.0f / 30;
    //if (inumiterSIMPLE371 > 30) {
      //  lower = free_debug_parametr1;
   // }
    //lower = 1.0f / 60; // достаточно одной тридцатой.
    //power_iters = 0;
   // scale = false;


    doublereal hi, lo;

    hi = spectral_radius(sl, slb, maxelm, maxbound);

    lo = hi * lower;


   

    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();







    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);

    // пороговое значение невязки
    doublereal eps = 1e-40;


    integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.


    doublereal* data_b = new doublereal[maxbound * 2];
    int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
    for (int i = 0; i < maxbound; ++i) {

        const int is = i * 2;
        data_b[is] = slb[i].ai;
        data_b[is + 1] = slb[i].aw;

        if (slb[i].iI > -1) {
            col_ind_b[is] = slb[i].iI;
        }
        else {
            data_b[is] = 0.0;
            col_ind_b[is] = slb[i].iW;
        }
        col_ind_b[is + 1] = slb[i].iW;
    }

    delete[] slb;

    // Улучшает попадание в кеш ускорение на 12.5%.
    doublereal* data = new doublereal[maxelm * 7];
    int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
    for (int i = 0; i < maxelm; ++i) {

        const int is = i * 7;
        data[is] = sl[i].ae;
        data[is + 1] = sl[i].aw;
        data[is + 2] = sl[i].an;
        data[is + 3] = sl[i].as;
        data[is + 4] = sl[i].at;
        data[is + 5] = sl[i].ab;
        data[is + 6] = sl[i].ap;

        col_ind[is] = sl[i].iE;
        col_ind[is + 1] = sl[i].iW;
        col_ind[is + 2] = sl[i].iN;
        col_ind[is + 3] = sl[i].iS;
        col_ind[is + 4] = sl[i].iT;
        col_ind[is + 5] = sl[i].iB;
        col_ind[is + 6] = sl[i].iP;
    }

    delete[] sl;


    const integer size0 = maxelm + maxbound;

    doublereal* x_new = new doublereal[size0];

    // Вектора необходимые для работы предобусловленного метода сопряжённых градиентов.
    doublereal* r75 = nullptr;
    doublereal* p75 = nullptr;
    doublereal* dx75 = nullptr;
    doublereal* dax75 = nullptr;
    doublereal* z76 = nullptr;
    doublereal* z77 = nullptr;
    doublereal* s76 = nullptr;
    doublereal* zold = nullptr;
    doublereal* znew = nullptr;

    r75 = new doublereal[size0];
    p75 = new doublereal[size0];
    dx75 = new doublereal[size0];
    dax75 = new doublereal[size0];
    z76 = new doublereal[size0];
    z77 = new doublereal[size0];
    s76 = new doublereal[size0];
    zold = new doublereal[size0];
    znew = new doublereal[size0];

#pragma omp parallel for
    for (integer i = 0; i < size0; ++i)
    {
        x_new[i] = 0.0;
       
        r75[i] = 0.0;
        p75[i] = 0.0;

        z76[i] = 0.0;
        s76[i] = 0.0;

        zold[i] = 0.0;
        znew[i] = 0.0;

        // результат умножения матрицы на вектор.
        dax75[i] = 0.0;
        // Начальное приближение.
        dx75[i] = x[i];
    }

    

    {
        // Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
        // опираемся на опыт теплопередачи.
        doublereal c1 = 0.026; // воздух
        doublereal c2 = 2000; // алмаз

        // см. Николаев
        if (bprintMessage) {
            std::cout << " number equations " << maxelm + maxbound << std::endl;
        }
        kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
        if (bprintMessage) {
            std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
        }
    }
    
   

    doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
    integer iout = 0;
 
   
    // результат записан в r75.
    // r75 - вектор невязки.
    residual_cg(
        maxelm, maxbound,
        data, col_ind,
        data_b, col_ind_b,
        dx75, b, r75);

#pragma omp parallel for 
    for (integer i75 = 0; i75 < size0; ++i75) {

        z76[i75] = 0.0;
        z77[i75] = 0.0;
        s76[i75] = r75[i75];
        // r75[i75] = s76[i75 + 1];
    }

    int iprecond = 0;

    if (iprecond==0)
    {
        bool bmem = false;

        doublereal r_start;

        const doublereal tau0 = 2.0 / (hi + lo);
        const doublereal ksi = lo / hi;
        const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


        //precontidioner
        // A, b, x -> x_new
        precontidioner(degree,
            tau0,
            rho0,
            mu,
            maxelm, maxbound,
            dmax, dmax2, delta, r_start,
            bprintMessage, bsolidT, j,
            iout, iVar, size0, 
            bmem, eps,
            data, col_ind,
            // A*z76=s76;
            data_b, col_ind_b, z77, s76, z76);

        

    }
    else {

        residual_cg(
            maxelm, maxbound,
            data, col_ind,
            data_b, col_ind_b,
            z77, s76, z76);
    }

#pragma omp parallel for 
    for (integer i75 = 0; i75 < size0; ++i75) {
        // Возвращаем результат.
        zold[i75] = z76[i75];
        p75[i75] = zold[i75];
    }

    doublereal alpha75, beta75;

    doublereal dnew = Scal(r75, zold, size0);

    // Используется для досрочного прерывания вычислительного процесса
       // как в алгоритме FGMRES Юсефа Саада и Мартина Г. Шульца.
    doublereal norma_b = NormaV_for_gmres(b, size0);

    

    iout = 0;

    bool bmem = true;
    doublereal dres0;

    do {

        


#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            zold[i75] = 0.0;
        }

        // нумерация  p75 и dax начинается с нуля.
        //Matrix_by_vector_for_Cheb  
        residual_cg(
            maxelm, maxbound,
            data, col_ind,
            data_b, col_ind_b,
            p75, zold, dax75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            dax75[i75] *= -1;// результат занесён в  dax75
        }


        doublereal scal_val1 = Scal(p75, dax75, size0);

        if ((fabs(dnew) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
            alpha75 = 1.0;
        }
        else if (fabs(dnew) < 1e-30) {
            alpha75 = 0.0;
        }
        else {
            alpha75 = dnew / scal_val1;
        }

#pragma omp parallel for
        for (integer i75 = 0; i75 < size0; ++i75)
        {
            dx75[i75] += alpha75 * p75[i75];
            r75[i75] -= alpha75 * dax75[i75];
        }


        doublereal deltai75 = NormaV_for_gmres(r75, size0);
        dres0 = deltai75 / norma_b;

        if (dres0 < 1.0e-7) break;

        if (iout > 2600* degree) break; // 15.02.2017

        
       



        // Предобуславливание.
#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75)
        {
            z76[i75] = 0.0;
            z77[i75] = 0.0;
            s76[i75] = r75[i75];
        }

        if (iprecond == 0) 
        {

            

            doublereal r_start;

            const doublereal tau0 = 2.0 / (hi + lo);
            const doublereal ksi = lo / hi;
            const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


            //precontidioner
            // A, b, x -> x_new
            // A*z76=s76;
            precontidioner(degree,
                tau0,
                rho0,
                mu,
                maxelm, maxbound,
                dmax, dmax2, delta, r_start,
                bprintMessage, bsolidT, j,
                iout, iVar, size0, 
                bmem, eps,
                data, col_ind,
                data_b, col_ind_b, z77, s76, z76);

           

            delta /= r_start;

            // Корректировка нижней границы спектра lo линейного оператора А.
            // Жуков, Феодоритова. 13.02.2022
            correct_lo(lo, bmem, degree, hi, delta, bprintMessage);


            //std::cout << "lo = " << lo << std::endl;
        }
        else {

            residual_cg(
                maxelm, maxbound,
                data, col_ind,
                data_b, col_ind_b,
                z77, s76, z76);
        }

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {

            // Возвращаем результат.
            znew[i75] = z76[i75];
        }

        doublereal dold = dnew;
        dnew = Scal(znew, r75, size0);

        //beta75= scr75 / scr75_old;
        //beta75 = Scal(r75, znew,n75) / Scal(rold75, zold, n75);

        if ((fabs(dnew) < 1e-30) && (fabs(dold) < 1e-30)) {
            beta75 = 1.0;
        }
        else if (fabs(dnew) < 1e-30) {
            beta75 = 0.0;
        }
        else {
            beta75 = dnew / dold;
        }



#pragma omp parallel for
        for (integer i75 = 0; i75 < size0; ++i75) {
            // p75[i75] = r75[i75] + beta75 * p75[i75];
            p75[i75] = znew[i75] + beta75 * p75[i75];
        }

        //dres0 = NormaV(r75, size0);
        dres0 = NormaV_for_gmres(r75, size0) / norma_b;

        std::cout << iout << "  " << dres0 << std::endl;



        //getchar();

        ++iout;
    } while (dres0 >= 1.0e-7);


#pragma omp parallel for
    for (integer i75 = 0; i75 < size0; ++i75) {
        x[i75] = dx75[i75];
    }

    for (int i = 0; i < maxelm; ++i) {
        const int is = i * 7;
        const int iP = col_ind[is + 6];
        x[iP] *= D[iP];
    }

#pragma omp parallel for 
    for (int i = 0; i < maxbound; ++i) {
        const int is = i * 2;
        const int iW = col_ind_b[is + 1];
        x[iW] *= D[iW] ;
    }

    delete[] data;
    delete[] col_ind;

    delete[] data_b;
    delete[] col_ind_b;

    std::cout << " " << lo << " " << j << " ";

    if (j > kend * degree - 100) {
        std::cout << "\n chebyshev divergence " << dmax << " ";
        system("pause");
    }

    delete[] x_new;
    // Освобождение оперативной памяти.
    delete[] r75;
    delete[] p75;
    delete[] dx75;
    delete[] dax75;
    delete[] znew;
    delete[] zold;

    delete[] z76;
    delete[] s76;
   

    // Реанимируем то что было выделено до решения СЛАУ,
   // а в момент решения СЛАУ мы экономили оперативную память.
    slb = new equation3D_bon[maxbound];
    sl = new equation3D[maxelm];

    if (bprintMessage) {
        if ((iVar == PAM) || (iVar == TEMP)) {
            printf("calc complete...\n");
            system("pause");
        }
    }
    //printf("4000 %e \n", dmax);
    //system("pause");

} // cg_chebyshev_3Dnow



// . 
// Сделан на структурированной сетке. 
// 05.02.2022.
// 
// Метод BiCGStab предобусловленный методом Чебышева. 12.02.2022.
bool bicgstab_chebyshev_3Dnow_cache(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
    doublereal*& x, int maxelm, int maxbound, integer iVar,
    bool bprintMessage, bool bsolidT) {
    // К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
    // скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
    //printf("chebyshev3Dnow incomming...\n"); // debug.
    //system("pause");

    // Параметр релаксации лежит в интервале от 0.0 до 2.0.
    // Верхняя релаксация способна существенно ускорить вычислительный процесс. 
    // Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
    // rURF=1.5 можно использовать в уравнении на поправку давления.
    // В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
    // к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
    // в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
    // можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
    // коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
    // поэтому наверное лучше не вводить коээффициент 
    // верхней релаксации равный 1.5 в уравнениях на скорость. 

     // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

    std::cout << "\nbicgstab_chebyshev_3Dnow_cache\n";
    // Надо запомнить и заморозить матрицу СЛАУ.

    equation3D* sl_copy = new equation3D[maxelm];
#pragma omp parallel for
    for (int i = 0; i < maxelm; ++i) {
        sl_copy[i] = sl[i];
    }
    equation3D_bon* slb_copy = new equation3D_bon[maxbound];
#pragma omp parallel for
    for (int i = 0; i < maxbound; ++i) {
        slb_copy[i] = slb[i];
    }


#ifdef _OPENMP 
    omp_set_num_threads(number_cores()); // установка числа потоков
#endif

   
   

    const bool bactive_ren = true;

    if (!b_setup_CathilMC_Temp) {
        // оптимизация попадания в кэш.
        renumerate_setup(sl, slb, maxelm, maxbound, new_number_CathilMC_Temp, new_number_internal_CathilMC_Temp, new_number_bound_CathilMC_Temp, rev_number_CathilMC_Temp, bactive_ren);
        b_setup_CathilMC_Temp = true;
    }

    renumerate_direct(sl, slb, maxelm, maxbound, x, b, new_number_CathilMC_Temp, new_number_internal_CathilMC_Temp, new_number_bound_CathilMC_Temp, rev_number_CathilMC_Temp,  bactive_ren);


    // Чистый cg лучше сходится без деления на диагональ.
    if (1 && bsolidT) {
        // диагональное предобуславливание
        // для температуры в твёрдом теле.
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae /= sl[i].ap;
            sl[i].aw /= sl[i].ap;
            sl[i].an /= sl[i].ap;
            sl[i].as /= sl[i].ap;
            sl[i].at /= sl[i].ap;
            sl[i].ab /= sl[i].ap;

            b[sl[i].iP] /= sl[i].ap;
            sl[i].ap = 1.0;
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            slb[i].ai /= slb[i].aw;
            b[slb[i].iW] /= slb[i].aw;

            slb[i].aw = 1.0;
        }
    }

    std::vector<double> D(static_cast<integer>(maxelm) + static_cast<integer>(maxbound), 1.0);

    if (0 && bsolidT) {

        // Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            D[sl[i].iP] = 1 / sqrt(fabs(sl[i].ap));
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            D[slb[i].iW] = 1 / sqrt(fabs(slb[i].aw));
        }

        // Then, apply the scaling in-place:
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae *= D[sl[i].iP] * D[sl[i].iE];
            sl[i].aw *= D[sl[i].iP] * D[sl[i].iW];
            sl[i].an *= D[sl[i].iP] * D[sl[i].iN];
            sl[i].as *= D[sl[i].iP] * D[sl[i].iS];
            sl[i].at *= D[sl[i].iP] * D[sl[i].iT];
            sl[i].ab *= D[sl[i].iP] * D[sl[i].iB];
            sl[i].ap *= D[sl[i].iP] * D[sl[i].iP];

            b[sl[i].iP] *= D[sl[i].iP];
        }
        
#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            if (slb[i].iI > -1) {
                slb[i].ai *= D[slb[i].iW] * D[slb[i].iI];
            }
            slb[i].aw *= D[slb[i].iW] * D[slb[i].iW];

            b[slb[i].iW] *= D[slb[i].iW];
        }



    }

    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
   // bool scale;


    higher = 1.0f;
    lower = 1.0f / 3.0f;// 1.0f / 30;
    //if (inumiterSIMPLE371 > 30) {
      //  lower = free_debug_parametr1;
   // }
    //lower = 1.0f / 60; // достаточно одной тридцатой.
    //power_iters = 0;
   // scale = false;


    doublereal hi, lo; // верхняя hi и нижняя lo границы спектра.

    // Вычисление верхней границы спектра по теореме Гершгорина о кругах.
    hi = spectral_radius(sl, slb, maxelm, maxbound);

    lo = hi * lower;


  

    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();







    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);

    // пороговое значение невязки
    doublereal eps = 1e-40;


    integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.

    if (iVar == TEMP) {
        if (maxbound < maxelm) {
            // Т.к. данный метод предназначен только для структурированной сетки то мы 
            // переливаем данные из slb на свободные места из-за АЛИС в sl не тратя ни грамма лишненй памяти.
#pragma omp parallel for 
            for (int i = 0; i < maxbound; ++i) {

                sl[i].ab2 = slb[i].ai;
                sl[i].ab3 = slb[i].aw;
                sl[i].iB2 = slb[i].iI;
                sl[i].iB3 = slb[i].iW;
            }

            // Сразу освобождаем оперативную память.
            delete[] slb;           

        }
    }

    doublereal* data_b = new doublereal[static_cast<integer>(maxbound) * 2];
    int* col_ind_b = new int[static_cast<integer>(maxbound) * 2];

    if ((iVar == PAM) || (maxbound >= maxelm)) {
#pragma omp parallel for
        for (int i = 0; i < maxbound; ++i) {

            const int is = i * 2;
            data_b[is] = slb[i].ai;
            data_b[is + 1] = slb[i].aw;

            if (slb[i].iI > -1) {
                col_ind_b[is] = slb[i].iI;
            }
            else {
                data_b[is] = 0.0;
                col_ind_b[is] = slb[i].iW;
            }
            col_ind_b[is + 1] = slb[i].iW;
        }

        // Сразу освобождаем оперативную память.
        if (iVar == TEMP) {
            delete[] slb;
        }
    }
    else {

#pragma omp parallel for
        for (int i = 0; i < maxbound; ++i) {

            const int is = i * 2;
            data_b[is] = sl[i].ab2;
            data_b[is + 1] = sl[i].ab3;

            if (sl[i].iB2 > -1) {
                col_ind_b[is] = sl[i].iB2;
            }
            else {
                data_b[is] = 0.0;
                col_ind_b[is] = sl[i].iB3;
            }
            col_ind_b[is + 1] = sl[i].iB3;
        }
    }


    // Улучшает попадание в кеш ускорение на 12.5%.
    doublereal* data = new doublereal[static_cast<integer>(maxelm) * 7];
    int* col_ind = new int[static_cast<integer>(maxelm) * 7];

#pragma omp parallel for 
    for (int i = 0; i < maxelm; ++i) {

        const int is = i * 7;
        data[is] = sl[i].ae;
        data[is + 1] = sl[i].aw;
        data[is + 2] = sl[i].an;
        data[is + 3] = sl[i].as;
        data[is + 4] = sl[i].at;
        data[is + 5] = sl[i].ab;
        data[is + 6] = sl[i].ap;

        col_ind[is] = sl[i].iE;
        col_ind[is + 1] = sl[i].iW;
        col_ind[is + 2] = sl[i].iN;
        col_ind[is + 3] = sl[i].iS;
        col_ind[is + 4] = sl[i].iT;
        col_ind[is + 5] = sl[i].iB;
        col_ind[is + 6] = sl[i].iP;
    }

    if (iVar == TEMP) {
        delete[] sl;
    }



    const integer size0 = static_cast<integer>(maxelm) + static_cast<integer>(maxbound);

    // 16 векторов это 2.3 матрицы СЛАУ.
    // 7 векторов дополнительная матрица для кеша.
    // Всего 3.3 матрицы. 23 вектора.


    // Вектора необходимые для работы предобусловленного метода сопряжённых градиентов.
   // Вектора необходимые для работы BiCGStab.
    doublereal* ri75 = nullptr;
    doublereal* roc75 = nullptr;
    doublereal* s75 = nullptr;
    doublereal* t75 = nullptr;
    doublereal* vec75 = nullptr;
    doublereal* vi75 = nullptr;
    doublereal* pi75 = nullptr;
    doublereal* dx75 = nullptr;

    doublereal* y75 = nullptr;
    doublereal* z75 = nullptr;
    // Первое предобуславливание:
    doublereal* pi76 = nullptr;

    // Второе предобуславливание:
    doublereal* z76 = nullptr;
    doublereal* s76 = nullptr;
    doublereal* z77 = nullptr;

    pi76 = new doublereal[size0];
    ri75 = new doublereal[size0];
    roc75 = new doublereal[size0];
    s75 = new doublereal[size0];
    t75 = new doublereal[size0];
    vec75 = new doublereal[size0];
    vi75 = new doublereal[size0];
    pi75 = new doublereal[size0];
    dx75 = new doublereal[size0];
    y75 = new doublereal[size0];
    z75 = new doublereal[size0];
    z77 = new doublereal[size0];
    z76 = new doublereal[size0];
    s76 = new doublereal[size0];
    

#pragma omp parallel for
    for (integer i = 0; i < size0; ++i)
    {

        s75[i] = 0.0;
        t75[i] = 0.0;

        vi75[i] = 0.0;
        pi75[i] = 0.0;

        y75[i] = 0.0;
        z75[i] = 0.0;

        // Начальное приближение.
        dx75[i] = x[i];
    }

    
    doublereal delta075 = 1.0e30;// deltai75 = 1.0e30;
    doublereal bet75 = 0.0, roi75 = 0.0;
    doublereal roim175 = 1.0, al75 = 1.0, wi75 = 1.0;

    {
        // Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
        // опираемся на опыт теплопередачи.
        doublereal c1 = 0.026; // воздух
        doublereal c2 = 2000; // алмаз

        // см. Николаев
        if (bprintMessage) {
            std::cout << " number equations " << maxelm + maxbound << std::endl;
        }
        kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
        if (bprintMessage) {
            std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
        }
    }





    doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
    integer iout = 0;


    // результат записан в r75.
    // r75 - вектор невязки.
    residual_cg(
        maxelm, maxbound,
        data, col_ind,
        data_b, col_ind_b,
        x, b, ri75);


#pragma omp parallel for 
    for (integer i75 = 0; i75 < size0; ++i75) {

        roc75[i75] = 1.0;

        s76[i75] = ri75[i75];
        // r75[i75] = s76[i75 + 1];
    }

    delta075 = NormaV(ri75, size0);
    if (iVar == TEMP) {
        std::wcout << delta075 << std::endl;
    }

    // Используется для досрочного прерывания вычислительного процесса
        // как в алгоритме FGMRES Юсефа Саада и Мартина Г. Шульца.
    doublereal norma_b = NormaV_for_gmres(b, size0);

    doublereal dres00 = (NormaV_for_gmres(ri75, size0) / norma_b);

    integer iprecond = 0;


    

    iout = 0;

    bool bmem = true;
    doublereal dres0;

    bool bdivergence_chebyshev = false;

    do {


        roi75 = Scal(roc75, ri75, size0);
        bet75 = (roi75 / roim175) * (al75 / wi75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75) * bet75;
            pi75[i75] = pibuf75;
        }


        // Первое предобуславливание.
        // Ky=pi
#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

           //y75[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7*(rand()%1000));


            z77[i75] = 0.0;
            pi76[i75] = pi75[i75];
        }

        if (iprecond == 0)
        {

            doublereal r_start;

            // do {



            const doublereal tau0 = 2.0 / (hi + lo);
            const doublereal ksi = lo / hi;
            const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


            //precontidioner
            // A, b, x -> z77
            // A*z77=pi76;
            precontidioner(degree,
                tau0,
                rho0,
                mu,
                maxelm, maxbound,
                dmax, dmax2, delta, r_start,
                bprintMessage, bsolidT, j,
                iout, iVar, size0,
                bmem, eps,
                data, col_ind,
                data_b, col_ind_b, y75, pi76, z77);



            delta /= r_start;

            // Корректировка нижней границы спектра lo линейного оператора А.
            // Жуков, Феодоритова. 13.02.2022
            correct_lo(lo, bmem, degree, hi, delta, bprintMessage);

            // std::wcout << dmax << std::endl;

         //} while (dmax >= 0.0001);
         //std::cout << "lo = " << lo << std::endl;

         //getchar();
        }
        else {

            residual_cg(
                maxelm, maxbound,
                data, col_ind,
                data_b, col_ind_b,
                y75, pi76, z77);
        }

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {

            // Возвращаем результат.            
            y75[i75] = z77[i75];
            z77[i75] = 0.0;
        }



        // vi==A*y;
        // нумерация  y75 и vi75 начинается с нуля.
        //Matrix_by_vector_for_Cheb  
        residual_cg(
            maxelm, maxbound,
            data, col_ind,
            data_b, col_ind_b,
            y75, z77, vi75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            vi75[i75] *= -1;// результат занесён в  dax75
        }


        doublereal scal_val1 = Scal(roc75, vi75, size0);

        if ((fabs(roi75) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
            al75 = 1.0;
        }
        else if (fabs(roi75) < 1e-30) {
            al75 = 0.0;
        }
        else {
            al75 = roi75 / scal_val1;
        }

#pragma omp parallel for
        for (integer i75 = 0; i75 < size0; ++i75)
        {
            s75[i75] = ri75[i75] - al75 * vi75[i75];
            z75[i75] = 0.0;
        }

        // Предобуславливание.
#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75)
        {
            z77[i75] = 0.0;// только ноль
            //z77[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7 * (rand() % 1000));


            vec75[i75] = s75[i75];
            z76[i75] = 0.0;
            s76[i75] = s75[i75];
        }

        // A*z76=s76;

        if (iprecond == 0)
        {

            doublereal r_start;

            // do {



            const doublereal tau0 = 2.0 / (hi + lo);
            const doublereal ksi = lo / hi;
            const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


            //std::wcout << Scal(s76, s76, size0) << std::endl;
           // getchar();

            //precontidioner
            // A, b, x -> z76
            // A*z76=s76;
            precontidioner(degree,
                tau0,
                rho0,
                mu,
                maxelm, maxbound,
                dmax, dmax2, delta, r_start,
                bprintMessage, bsolidT, j,
                iout, iVar, size0,
                bmem, eps,
                data, col_ind,
                data_b, col_ind_b, z77, s76, z76);



            delta /= r_start;

            // Корректировка нижней границы спектра lo линейного оператора А.
            // Жуков, Феодоритова. 13.02.2022
            correct_lo(lo, bmem, degree, hi, delta, bprintMessage);


            // std::wcout << dmax << std::endl;

            // } while (dmax >= 0.0001);

             //getchar();

             //std::cout << "lo = " << lo << std::endl;
        }
        else {

            residual_cg(
                maxelm, maxbound,
                data, col_ind,
                data_b, col_ind_b,
                z77, s76, z76);
        }

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {

            z77[i75] = 0.0;
            // Возвращаем результат.
            s75[i75] = vec75[i75];
            // Возвращаем результат.
            z75[i75] = z76[i75];
        }

        // t==A*z;
         // нумерация  y75 и vi75 начинается с нуля.
        //Matrix_by_vector_for_Cheb  
        residual_cg(
            maxelm, maxbound,
            data, col_ind,
            data_b, col_ind_b,
            z75, z77, t75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            t75[i75] *= -1;// результат занесён в  dax75
        }


        wi75 = Scal(t75, s75, size0) / Scal(t75, t75, size0);
        // std::cout << "Scal(t75,s75,n75)==" << Scal(t75,s75,n75) << ", Scal(t75,t75,n75)=="<< Scal(t75,t75,n75) << std::endl;

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            //dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
            dx75[i75] += al75 * y75[i75] + wi75 * z75[i75]; // так стало с предобуславливателем
            ri75[i75] = s75[i75] - wi75 * t75[i75];
        }
        //dres0 = deltai75 = NormaV(ri75, size0);

        //dres0 = (NormaV_for_gmres(ri75, size0) / norma_b)/ size0;
        dres0 = (NormaV_for_gmres(ri75, size0) / norma_b);
        //dres0 = (NormaV_for_gmres(ri75, size0));

        if (iVar == TEMP) {
            std::cout << iout << "  " << dres0 << std::endl;
        }

        if (iVar == PAM) {
            // CFD break for PAM.
            if (((dres0 / dres00 < 0.01) && ((size0 < 6000000) || (iout > 3)))) break;

            if (iout > (size0 > 6000000 ? 15 : 5)) {
                if (dres0 / dres00 < (size0 > 6000000 ? 0.03 : 0.1)) {
                    break;
                }
                else {
                    if (iout > (size0 > 6000000 ? 21 : 7)) {
                        std::cout << " PAM chebyshev degree = " << degree << " divergence!#it = " << (size0 > 6000000 ? 21 : 7) << "\n";
                        //system("PAUSE");
                        bdivergence_chebyshev = true;
                        break;
                    }
                }
            }
        }

        //getchar();

        ++iout;
    } while (dres0 >= 1.0e-7);


    if (iVar == PAM) {
        std::cout << " " << iout << " ";
    }

    if (!bdivergence_chebyshev) {

#pragma omp parallel for
        for (integer i75 = 0; i75 < size0; ++i75) {

            x[i75] = dx75[i75];
        }



        for (int i = 0; i < maxelm; ++i) {
            const int is = i * 7;
            const int iP = col_ind[is + 6];
            x[iP] *= D[iP];
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            const int is = i * 2;
            const int iW = col_ind_b[is + 1];
            x[iW] *= D[iW];
        }
    }

    delete[] data;
    delete[] col_ind;

    delete[] data_b;
    delete[] col_ind_b;

    if (iVar == TEMP) {
        std::cout << " " << lo << " " << j << " ";
    }

    if (j > kend * degree - 100) {
        std::cout << "\n chebyshev divergence " << dmax << " ";
        system("pause");
    }


    // Освобождение оперативной памяти.
    // Первое предобуславливание

    delete[] pi76;
    pi76 = nullptr;


    // Второе предобуславливание

    delete[] z76;
    z76 = nullptr;
    delete[] s76;
    s76 = nullptr;

    delete[] ri75;
    ri75 = nullptr;
    delete[] roc75;
    roc75 = nullptr;

    delete[] s75;
    s75 = nullptr;
    delete[] t75;
    t75 = nullptr;

    delete[] vec75;
    vec75 = nullptr;
    delete[] vi75;
    vi75 = nullptr;


    delete[] pi75;
    pi75 = nullptr;
    delete[] dx75;
    dx75 = nullptr;


    delete[] y75;
    y75 = nullptr;

    delete[] z75;
    z75 = nullptr;

    delete[] z77;
    z77 = nullptr;

    if (iVar == TEMP) {
        // Реанимируем то что было выделено до решения СЛАУ,
        // а в момент решения СЛАУ мы экономили оперативную память.
        slb = new equation3D_bon[maxbound];
        sl = new equation3D[maxelm];
    }

    // Восстанавливаем исходную нумерацию.
    renumerate_reverse(x, rev_number_CathilMC_Temp, maxelm, maxbound, bactive_ren);
    renumerate_reverse(b, rev_number_CathilMC_Temp, maxelm, maxbound, bactive_ren);

#pragma omp parallel for
    for (int i = 0; i < maxelm; ++i) {
        sl[i] = sl_copy[i];
    }
  
    delete[] sl_copy;

#pragma omp parallel for
    for (int i = 0; i < maxbound; ++i) {
         slb[i] = slb_copy[i];
    }

    delete[] slb_copy;


    // Восстанавливаем исходную нумерацию т.к. метод Чебышева мог несойтись и для последующего старта amg нам нужны нетронутые данные.


    if (bprintMessage) {
        if ((iVar == PAM) || (iVar == TEMP)) {
            printf("calc complete...\n");
            system("pause");
        }
    }
    //printf("4000 %e \n", dmax);
    //system("pause");

    return bdivergence_chebyshev;

} // bicgstab_chebyshev_3Dnow_cache

// . 
// Сделан на структурированной сетке. 
// 05.02.2022.
// 
// Метод BiCGStab предобусловленный методом Чебышева. 12.02.2022.
bool bicgstab_chebyshev_3Dnow(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
    doublereal*& x, integer maxelm, integer maxbound, integer iVar,
    bool bprintMessage, bool bsolidT) {
    // К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
    // скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
    //printf("chebyshev3Dnow incomming...\n"); // debug.
    //system("pause");

    // Параметр релаксации лежит в интервале от 0.0 до 2.0.
    // Верхняя релаксация способна существенно ускорить вычислительный процесс. 
    // Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
    // rURF=1.5 можно использовать в уравнении на поправку давления.
    // В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
    // к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
    // в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
    // можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
    // коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
    // поэтому наверное лучше не вводить коээффициент 
    // верхней релаксации равный 1.5 в уравнениях на скорость. 

     // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

#ifdef _OPENMP 
    omp_set_num_threads(number_cores()); // установка числа потоков
#endif

    
    // Чистый cg лучше сходится без деления на диагональ.
    if (1&&bsolidT) {
        // диагональное предобуславливание
        // для температуры в твёрдом теле.
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae /= sl[i].ap;
            sl[i].aw /= sl[i].ap;
            sl[i].an /= sl[i].ap;
            sl[i].as /= sl[i].ap;
            sl[i].at /= sl[i].ap;
            sl[i].ab /= sl[i].ap;

            b[sl[i].iP] /= sl[i].ap;
            sl[i].ap = 1.0;
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            slb[i].ai /= slb[i].aw;
            b[slb[i].iW] /= slb[i].aw;

            slb[i].aw = 1.0;
        }
    }

    std::vector<double> D(maxelm + maxbound, 1.0);

    if (0&&bsolidT) {

        // Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            D[sl[i].iP] = 1 / sqrt(fabs(sl[i].ap));
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            D[slb[i].iW] = 1 / sqrt(fabs(slb[i].aw));
        }

        // Then, apply the scaling in-place:
#pragma omp parallel for 
        for (int i = 0; i < maxelm; ++i) {
            sl[i].ae *= D[sl[i].iP] * D[sl[i].iE];
            sl[i].aw *= D[sl[i].iP] * D[sl[i].iW];
            sl[i].an *= D[sl[i].iP] * D[sl[i].iN];
            sl[i].as *= D[sl[i].iP] * D[sl[i].iS];
            sl[i].at *= D[sl[i].iP] * D[sl[i].iT];
            sl[i].ab *= D[sl[i].iP] * D[sl[i].iB];
            sl[i].ap *= D[sl[i].iP] * D[sl[i].iP];

            b[sl[i].iP] *= D[sl[i].iP];
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            if (slb[i].iI > -1) {
                slb[i].ai *= D[slb[i].iW] * D[slb[i].iI];
            }
            slb[i].aw *= D[slb[i].iW] * D[slb[i].iW];

            b[slb[i].iW] *= D[slb[i].iW];
        }



    }
    
    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
   // bool scale;


    higher = 1.0f;
    lower = 1.0f / 3.0f;// 1.0f / 30;
    //if (inumiterSIMPLE371 > 30) {
      //  lower = free_debug_parametr1;
   // }
    //lower = 1.0f / 60; // достаточно одной тридцатой.
    //power_iters = 0;
   // scale = false;


    doublereal hi, lo; // верхняя hi и нижняя lo границы спектра.

    // Вычисление верхней границы спектра по теореме Гершгорина о кругах.
    hi = spectral_radius(sl, slb, maxelm, maxbound);

    lo = hi * lower;


    //doublereal rURF = 1.0; // параметр верхней релаксации

    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();







    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);

    // пороговое значение невязки
    doublereal eps = 1e-40;


    integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.

    if (iVar == TEMP) {
        if (maxbound < maxelm) {
            // Т.к. данный метод предназначен только для структурированной сетки то мы 
            // переливаем данные из slb на свободные места из-за АЛИС в sl не тратя ни грамма лишненй памяти.
#pragma omp parallel for 
            for (int i = 0; i < maxbound; ++i) {

                sl[i].ab2 = slb[i].ai;
                sl[i].ab3 = slb[i].aw;
                sl[i].iB2 = slb[i].iI;
                sl[i].iB3 = slb[i].iW;
            }

            // Сразу освобождаем оперативную память.
            delete[] slb;            

        }
    }

    doublereal* data_b = new doublereal[maxbound * 2];
    int* col_ind_b = new int[maxbound * 2];

    if ((iVar==PAM) || (maxbound >= maxelm)) {
       #pragma omp parallel for
            for (int i = 0; i < maxbound; ++i) {

                const int is = i * 2;
                data_b[is] = slb[i].ai;
                data_b[is + 1] = slb[i].aw;

                if (slb[i].iI > -1) {
                    col_ind_b[is] = slb[i].iI;
                }
                else {
                    data_b[is] = 0.0;
                    col_ind_b[is] = slb[i].iW;
                }
                col_ind_b[is + 1] = slb[i].iW;
            }

            // Сразу освобождаем оперативную память.
            if (iVar == TEMP) {
                delete[] slb;
            }
    }
    else {

#pragma omp parallel for
        for (int i = 0; i < maxbound; ++i) {

            const int is = i * 2;
            data_b[is] = sl[i].ab2;
            data_b[is + 1] = sl[i].ab3;

            if (sl[i].iB2 > -1) {
                col_ind_b[is] = sl[i].iB2;
            }
            else {
                data_b[is] = 0.0;
                col_ind_b[is] = sl[i].iB3;
            }
            col_ind_b[is + 1] = sl[i].iB3;
        }
    }
    

     // Улучшает попадание в кеш ускорение на 12.5%.
    doublereal* data = new doublereal[maxelm * 7];
    int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
    for (int i = 0; i < maxelm; ++i) {

        const int is = i * 7;
        data[is] = sl[i].ae;
        data[is + 1] = sl[i].aw;
        data[is + 2] = sl[i].an;
        data[is + 3] = sl[i].as;
        data[is + 4] = sl[i].at;
        data[is + 5] = sl[i].ab;
        data[is + 6] = sl[i].ap;

        col_ind[is] = sl[i].iE;
        col_ind[is + 1] = sl[i].iW;
        col_ind[is + 2] = sl[i].iN;
        col_ind[is + 3] = sl[i].iS;
        col_ind[is + 4] = sl[i].iT;
        col_ind[is + 5] = sl[i].iB;
        col_ind[is + 6] = sl[i].iP;
    }

    if (iVar == TEMP) {
        delete[] sl;
    }
    


    const integer size0 = maxelm + maxbound;

    // 16 векторов это 2.3 матрицы СЛАУ.
    // 7 векторов дополнительная матрица для кеша.
    // Всего 3.3 матрицы. 23 вектора.
    

    // Вектора необходимые для работы предобусловленного метода сопряжённых градиентов.
   // Вектора необходимые для работы BiCGStab.
    doublereal* ri75 = nullptr;
    doublereal* roc75 = nullptr;
    doublereal* s75 = nullptr;
    doublereal* t75 = nullptr;
    doublereal* vec75 = nullptr;
    doublereal* vi75 = nullptr;
    doublereal* pi75 = nullptr;
    doublereal* dx75 = nullptr;
   
    doublereal* y75 = nullptr;
    doublereal* z75 = nullptr;
    // Первое предобуславливание:
    doublereal* pi76 = nullptr;
    
    // Второе предобуславливание:
    doublereal* z76 = nullptr;
    doublereal* s76 = nullptr;    
    doublereal* z77 = nullptr;    

    pi76 = new doublereal[size0];
    ri75 = new doublereal[size0];
    roc75 = new doublereal[size0];
    s75 = new doublereal[size0];
    t75 = new doublereal[size0];
    vec75 = new doublereal[size0];
    vi75 = new doublereal[size0];
    pi75 = new doublereal[size0];
    dx75 = new doublereal[size0];   
    y75 = new doublereal[size0];
    z75 = new doublereal[size0];
    z77 = new doublereal[size0];
    z76 = new doublereal[size0];
    s76 = new doublereal[size0];
    

#pragma omp parallel for
    for (integer i = 0; i < size0; ++i)
    {       

        s75[i] = 0.0;
        t75[i] = 0.0;

        vi75[i] = 0.0;
        pi75[i] = 0.0;

        y75[i] = 0.0;
        z75[i] = 0.0;
        
        // Начальное приближение.
        dx75[i] = x[i];
    }

    
    doublereal delta075 = 1.0e30;// deltai75 = 1.0e30;
    doublereal bet75 = 0.0, roi75 = 0.0;
    doublereal roim175 = 1.0, al75 = 1.0, wi75 = 1.0;

    {
        // Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
        // опираемся на опыт теплопередачи.
        doublereal c1 = 0.026; // воздух
        doublereal c2 = 2000; // алмаз

        // см. Николаев
        if (bprintMessage) {
            std::cout << " number equations " << maxelm + maxbound << std::endl;
        }
        kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
        if (bprintMessage) {
            std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
        }
    }



   

    doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
    integer iout = 0;


    // результат записан в r75.
    // r75 - вектор невязки.
    residual_cg(
        maxelm, maxbound,
        data, col_ind,
        data_b, col_ind_b,
        x, b, ri75);

#pragma omp parallel for 
    for (integer i75 = 0; i75 < size0; ++i75) {

        roc75[i75] = 1.0;

        s76[i75] = ri75[i75];
        // r75[i75] = s76[i75 + 1];
    }

    delta075 = NormaV(ri75, size0);
    if (iVar == TEMP) {
        std::wcout << delta075 << std::endl;
    }

    // Используется для досрочного прерывания вычислительного процесса
        // как в алгоритме FGMRES Юсефа Саада и Мартина Г. Шульца.
    doublereal norma_b = NormaV_for_gmres(b, size0);

    doublereal dres00 = (NormaV_for_gmres(ri75, size0) / norma_b);

    integer iprecond = 0;


    

    iout = 0;

    bool bmem = true;
    doublereal dres0;

    bool bdivergence_chebyshev = false;

    do {
        

        roi75 = Scal(roc75, ri75, size0);
        bet75 = (roi75 / roim175) * (al75 / wi75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75) * bet75;
            pi75[i75] = pibuf75;
        }


        // Первое предобуславливание.
        // Ky=pi
#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
             y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
            
            //y75[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7*(rand()%1000));


            z77[i75] = 0.0;
            pi76[i75] = pi75[i75];
        }

        if (iprecond == 0)
        {

            doublereal r_start;

           // do {

                

                const doublereal tau0 = 2.0 / (hi + lo);
                const doublereal ksi = lo / hi;
                const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


                //precontidioner
                // A, b, x -> z77
                // A*z77=pi76;
                precontidioner(degree,
                    tau0,
                    rho0,
                    mu,
                    maxelm, maxbound,
                    dmax, dmax2, delta, r_start,
                    bprintMessage, bsolidT, j,
                    iout, iVar, size0, 
                    bmem, eps,
                    data, col_ind,
                    data_b, col_ind_b, y75, pi76, z77);



                delta /= r_start;

                // Корректировка нижней границы спектра lo линейного оператора А.
                // Жуков, Феодоритова. 13.02.2022
                correct_lo(lo, bmem, degree, hi, delta, bprintMessage);

               // std::wcout << dmax << std::endl;

            //} while (dmax >= 0.0001);
            //std::cout << "lo = " << lo << std::endl;

            //getchar();
        }
        else {

            residual_cg(
                maxelm, maxbound,
                data, col_ind,
                data_b, col_ind_b,
                y75, pi76, z77);
        }

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {

            // Возвращаем результат.            
            y75[i75] = z77[i75];
            z77[i75] = 0.0;
        }



        // vi==A*y;
        // нумерация  y75 и vi75 начинается с нуля.
        //Matrix_by_vector_for_Cheb  
        residual_cg(
            maxelm, maxbound,
            data, col_ind,
            data_b, col_ind_b,
            y75, z77, vi75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            vi75[i75] *= -1;// результат занесён в  dax75
        }


        doublereal scal_val1 = Scal(roc75, vi75, size0);

        if ((fabs(roi75) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
            al75 = 1.0;
        }
        else if (fabs(roi75) < 1e-30) {
            al75 = 0.0;
        }
        else {
            al75 = roi75 / scal_val1;
        }

#pragma omp parallel for
        for (integer i75 = 0; i75 < size0; ++i75)
        {
            s75[i75] = ri75[i75] - al75 * vi75[i75];
            z75[i75] = 0.0;            
        }

        // Предобуславливание.
#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75)
        {
            z77[i75] = 0.0;// только ноль
            //z77[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7 * (rand() % 1000));


            vec75[i75] = s75[i75];
            z76[i75] = 0.0;
            s76[i75] = s75[i75];
        }

        // A*z76=s76;

        if (iprecond == 0)
        {

            doublereal r_start;

           // do {

            

            const doublereal tau0 = 2.0 / (hi + lo);
            const doublereal ksi = lo / hi;
            const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


            //std::wcout << Scal(s76, s76, size0) << std::endl;
           // getchar();

            //precontidioner
            // A, b, x -> z76
            // A*z76=s76;
            precontidioner(degree,
                tau0,
                rho0,
                mu,
                maxelm, maxbound,
                dmax, dmax2, delta, r_start,
                bprintMessage, bsolidT, j,
                iout, iVar, size0, 
                bmem, eps,
                data, col_ind,
                data_b, col_ind_b, z77, s76, z76);



            delta /= r_start;

            // Корректировка нижней границы спектра lo линейного оператора А.
            // Жуков, Феодоритова. 13.02.2022
            correct_lo(lo, bmem,  degree, hi, delta,  bprintMessage);
           

           // std::wcout << dmax << std::endl;

           // } while (dmax >= 0.0001);

            //getchar();

            //std::cout << "lo = " << lo << std::endl;
        }
        else {

            residual_cg(
                maxelm, maxbound,
                data, col_ind,
                data_b, col_ind_b,
                z77, s76, z76);
        }

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {

            z77[i75] = 0.0;
            // Возвращаем результат.
            s75[i75] = vec75[i75];
            // Возвращаем результат.
            z75[i75] = z76[i75];
        }

        // t==A*z;
         // нумерация  y75 и vi75 начинается с нуля.
        //Matrix_by_vector_for_Cheb  
        residual_cg(
            maxelm, maxbound,
            data, col_ind,
            data_b, col_ind_b,
            z75, z77, t75);

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            t75[i75] *= -1;// результат занесён в  dax75
        }


        wi75 = Scal(t75, s75, size0) / Scal(t75, t75, size0);
        // std::cout << "Scal(t75,s75,n75)==" << Scal(t75,s75,n75) << ", Scal(t75,t75,n75)=="<< Scal(t75,t75,n75) << std::endl;

#pragma omp parallel for 
        for (integer i75 = 0; i75 < size0; ++i75) {
            //dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
            dx75[i75] += al75 * y75[i75] + wi75 * z75[i75]; // так стало с предобуславливателем
            ri75[i75] = s75[i75] - wi75 * t75[i75];
        }
        //dres0 = deltai75 = NormaV(ri75, size0);

        //dres0 = (NormaV_for_gmres(ri75, size0) / norma_b)/ size0;
        dres0 = (NormaV_for_gmres(ri75, size0) / norma_b);
        //dres0 = (NormaV_for_gmres(ri75, size0));

        if (iVar == TEMP) {
            std::cout << iout << "  " << dres0 << std::endl;
        }

        if (iVar == PAM) {
            // CFD break for PAM.
            if (((dres0 / dres00 < 0.01)&&((size0 < 6000000)||(iout>3)))) break;

            if (iout > (size0 > 6000000 ? 15 : 5)) {
                if (dres0 / dres00 < (size0 > 6000000 ? 0.03 : 0.1)) {
                    break;
                }
                else {
                    if (iout > (size0 > 6000000 ? 21 : 7)) {
                        std::cout << " PAM chebyshev degree = " << degree << " divergence!#it = " << (size0 > 6000000 ? 21 : 7) << "\n";
                        //system("PAUSE");
                        bdivergence_chebyshev = true;
                        break;
                    }
                }
            }
        }

        //getchar();

        ++iout;
    } while (dres0 >= 1.0e-7);


    if (iVar == PAM) {
        std::cout << " " << iout << " ";
    }

    if (!bdivergence_chebyshev) {

#pragma omp parallel for
        for (integer i75 = 0; i75 < size0; ++i75) {

            x[i75] = dx75[i75];
        }



        for (int i = 0; i < maxelm; ++i) {
            const int is = i * 7;
            const int iP = col_ind[is + 6];
            x[iP] *= D[iP];
        }

#pragma omp parallel for 
        for (int i = 0; i < maxbound; ++i) {
            const int is = i * 2;
            const int iW = col_ind_b[is + 1];
            x[iW] *= D[iW];
        }
    }

    delete[] data;
    delete[] col_ind;

    delete[] data_b;
    delete[] col_ind_b;

    if (iVar == TEMP) {
        std::cout << " " << lo << " " << j << " ";
    }

    if (j > kend * degree - 100) {
        std::cout << "\n chebyshev divergence " << dmax << " ";
        system("pause");
    }

    
    // Освобождение оперативной памяти.
    // Первое предобуславливание

    delete[] pi76;
    pi76 = nullptr;
    

    // Второе предобуславливание

    delete[] z76;
    z76 = nullptr;
    delete[] s76;
    s76 = nullptr;

    delete[] ri75;
    ri75 = nullptr;
    delete[] roc75;
    roc75 = nullptr;

    delete[] s75;
    s75 = nullptr;
    delete[] t75;
    t75 = nullptr;

    delete[] vec75;
    vec75 = nullptr;
    delete[] vi75;
    vi75 = nullptr;


    delete[] pi75;
    pi75 = nullptr;
    delete[] dx75;
    dx75 = nullptr;


    delete[] y75;
    y75 = nullptr;

    delete[] z75;
    z75 = nullptr;

    delete[] z77;
    z77 = nullptr;

    if (iVar == TEMP) {
        // Реанимируем то что было выделено до решения СЛАУ,
        // а в момент решения СЛАУ мы экономили оперативную память.
        slb = new equation3D_bon[maxbound];
        sl = new equation3D[maxelm];
    }

    if (bprintMessage) {
        if ((iVar == PAM) || (iVar == TEMP)) {
            printf("calc complete...\n");
            system("pause");
        }
    }
    //printf("4000 %e \n", dmax);
    //system("pause");

    return bdivergence_chebyshev;

} // bicgstab_chebyshev_3Dnow

// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.

doublereal spectral_radius_Alice(equation3D*& sl, equation3D_bon*& slb,
    integer maxelm, integer maxbound)
{
    // Теорема Гершгорина о кругах.

    doublereal radius;

    //if (power_iters <= 0)
    {
        // Use Gershgorin disk theorem.
        radius = 0;

#pragma omp parallel
        {
            doublereal emax = 0;
            // doublereal  dia = 1.0;



#pragma omp for nowait
            for (integer i = 0; i < maxelm; ++i) {

                doublereal s = 0;


                s += ((sl[i].iE > -1) ? sl[i].ae : 0.0) + ((sl[i].iW > -1) ? sl[i].aw : 0.0) + 
                    ((sl[i].iN > -1) ? sl[i].an : 0.0) + ((sl[i].iS > -1) ? sl[i].as : 0.0) + 
                    ((sl[i].iT > -1) ? sl[i].at : 0.0) + ((sl[i].iB > -1) ? sl[i].ab : 0.0) +  sl[i].ap;

                s += ((sl[i].iE2 > -1) ? sl[i].ae2 : 0.0) + ((sl[i].iW2 > -1) ? sl[i].aw2 : 0.0) +
                    ((sl[i].iN2 > -1) ? sl[i].an2 : 0.0) + ((sl[i].iS2 > -1) ? sl[i].as2 : 0.0) +
                    ((sl[i].iT2 > -1) ? sl[i].at2 : 0.0) + ((sl[i].iB2 > -1) ? sl[i].ab2 : 0.0);

                s += ((sl[i].iE3 > -1) ? sl[i].ae3 : 0.0) + ((sl[i].iW3 > -1) ? sl[i].aw3 : 0.0) +
                    ((sl[i].iN3 > -1) ? sl[i].an3 : 0.0) + ((sl[i].iS3 > -1) ? sl[i].as3 : 0.0) +
                    ((sl[i].iT3 > -1) ? sl[i].at3 : 0.0) + ((sl[i].iB3 > -1) ? sl[i].ab3 : 0.0);

                s += ((sl[i].iE4 > -1) ? sl[i].ae4 : 0.0) + ((sl[i].iW4 > -1) ? sl[i].aw4 : 0.0) +
                    ((sl[i].iN4 > -1) ? sl[i].an4 : 0.0) + ((sl[i].iS4 > -1) ? sl[i].as4 : 0.0) +
                    ((sl[i].iT4 > -1) ? sl[i].at4 : 0.0) + ((sl[i].iB4 > -1) ? sl[i].ab4 : 0.0);
               
                //dia = sl[i].ap;

                emax = std::max(emax, s);
            }

#pragma omp for nowait
            for (integer i = 0; i < maxbound; ++i) {

                doublereal s = 0;

                if (slb[i].iI == -1) {
                    // s += fabs(slb[i].aw);
                    s += slb[i].aw;
                }
                else {
                    //s += fabs(slb[i].aw) + fabs(slb[i].ai);
                    s += slb[i].aw + slb[i].ai;
                }
                //dia = slb[i].aw;

                emax = std::max(emax, s);
            }

#pragma omp critical
            radius = std::max(radius, emax);
        }
    }
    //std::cout << "spectral radius" << std::endl;

    return  (radius < 0 ? static_cast<doublereal>(2.0) : radius);
}


// применяется для уравнения поправки давления. 
// Работает только на АЛИС сетке данная версия.
// 07.02.2022.
// 
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
void chebyshev3Dnow_Alice(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
    doublereal*& x, integer maxelm, integer maxbound, integer iVar,
    bool bprintMessage) {
    // К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
    // скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
    //printf("chebyshev3Dnow incomming...\n"); // debug.
    //system("pause");

    // Параметр релаксации лежит в интервале от 0.0 до 2.0.
    // Верхняя релаксация способна существенно ускорить вычислительный процесс. 
    // Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
    // rURF=1.5 можно использовать в уравнении на поправку давления.
    // В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
    // к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
    // в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
    // можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
    // коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
    // поэтому наверное лучше не вводить коээффициент 
    // верхней релаксации равный 1.5 в уравнениях на скорость. 

     // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

#ifdef _OPENMP 
    omp_set_num_threads(number_cores()); // установка числа потоков
#endif

    /// Chebyshev polynomial degree.
    const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

    /// highest eigen value safety upscaling.
    // use boosting factor for a more conservative upper bound estimate
    // See: Adams, Brezina, Hu, Tuminaro,
    //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
    //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
    //
    float higher;

    /// Lowest-to-highest eigen value ratio.
    float lower;

    // Number of power iterations to apply for the spectral radius
    // estimation. When 0, use Gershgorin disk theorem to estimate
    // spectral radius.
    //integer power_iters;

    // Scale the system matrix
   // bool scale;


    higher = 1.0f;
    lower = 1.0f / 3.0f;// 1.0f / 30;
    //if (inumiterSIMPLE371 > 30) {
      //  lower = free_debug_parametr1;
   // }
    //lower = 1.0f / 60; // достаточно одной тридцатой.
    //power_iters = 0;
   // scale = false;


    doublereal hi, lo;

    hi = spectral_radius_Alice(sl, slb, maxelm, maxbound);

    lo = hi * lower;


    doublereal rURF = 1.0; // параметр верхней релаксации

    hi *= higher;

    // Centre of ellipse containing the eigenvalues of A:
    //doublereal d = 0.5 * (hi + lo);

    // Semi-major axis of ellipse containing the eigenvalues of A:
    //doublereal c = 0.5 * (hi - lo);


    //static const scalar_type one = math::identity<scalar_type>();
    //static const scalar_type zero = math::zero<scalar_type>();







    // Переупорядоченное множество корней многочлена П.Л. Чебышева
    doublereal* mu = new doublereal[degree];
    Lebedev_Samarskii_Nikolaev(degree, mu);

    // пороговое значение невязки
    doublereal eps = 1e-40;


    integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.


    const integer size0 = maxelm + maxbound;

    doublereal* x_new = new doublereal[size0];
#pragma omp parallel for
    for (integer i = 0; i < size0; ++i)
    {
        x_new[i] = 0.0;
    }

    {
        // Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
        // опираемся на опыт теплопередачи.
        doublereal c1 = 0.026; // воздух
        doublereal c2 = 2000; // алмаз

        // см. Николаев
        if (bprintMessage) {
            std::cout << " number equations " << maxelm + maxbound << std::endl;
        }
        kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
        if (bprintMessage) {
            std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
        }
    }

   


    doublereal* data_b = new doublereal[maxbound * 2];
    int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
    for (int i = 0; i < maxbound; ++i) {

        const int is = i * 2;
        data_b[is] = slb[i].ai;
        data_b[is + 1] = slb[i].aw;

        if (slb[i].iI > -1) {
            col_ind_b[is] = slb[i].iI;
        }
        else {
            data_b[is] = 0.0;
            col_ind_b[is] = slb[i].iW;
        }
        col_ind_b[is + 1] = slb[i].iW;
    }

    doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
    integer iout = 0;
    while ((dmax > eps) && (j < kend * degree)) {


        bool bmem = true;

        doublereal r_start;

        const doublereal tau0 = 2.0 / (hi + lo);
        const doublereal ksi = lo / hi;
        const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

        for (int k = 0; k < degree; ++k) {



            rURF = tau0 / (1.0 + rho0 * mu[k]);

            dmax = 0.0;
            dmax2 = 0.0;
            delta = 0.0;

            doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
            for (int i = 0; i < maxelm; ++i) {

                doublereal sE, sW, sN, sS, sT, sB;

                if (sl[i].iE > -1) sE = sl[i].ae * x[sl[i].iE]; else sE = 0.0;
                if (sl[i].iW > -1) sW = sl[i].aw * x[sl[i].iW]; else sW = 0.0;
                if (sl[i].iN > -1) sN = sl[i].an * x[sl[i].iN]; else sN = 0.0;
                if (sl[i].iS > -1) sS = sl[i].as * x[sl[i].iS]; else sS = 0.0;
                if (sl[i].iT > -1) sT = sl[i].at * x[sl[i].iT]; else sT = 0.0;
                if (sl[i].iB > -1) sB = sl[i].ab * x[sl[i].iB]; else sB = 0.0;

                doublereal sE2, sW2, sN2, sS2, sT2, sB2;

                if (sl[i].iE2 > -1) sE2 = sl[i].ae2 * x[sl[i].iE2]; else sE2 = 0.0;
                if (sl[i].iW2 > -1) sW2 = sl[i].aw2 * x[sl[i].iW2]; else sW2 = 0.0;
                if (sl[i].iN2 > -1) sN2 = sl[i].an2 * x[sl[i].iN2]; else sN2 = 0.0;
                if (sl[i].iS2 > -1) sS2 = sl[i].as2 * x[sl[i].iS2]; else sS2 = 0.0;
                if (sl[i].iT2 > -1) sT2 = sl[i].at2 * x[sl[i].iT2]; else sT2 = 0.0;
                if (sl[i].iB2 > -1) sB2 = sl[i].ab2 * x[sl[i].iB2]; else sB2 = 0.0;

                doublereal sE3, sW3, sN3, sS3, sT3, sB3;

                if (sl[i].iE3 > -1) sE3 = sl[i].ae3 * x[sl[i].iE3]; else sE3 = 0.0;
                if (sl[i].iW3 > -1) sW3 = sl[i].aw3 * x[sl[i].iW3]; else sW3 = 0.0;
                if (sl[i].iN3 > -1) sN3 = sl[i].an3 * x[sl[i].iN3]; else sN3 = 0.0;
                if (sl[i].iS3 > -1) sS3 = sl[i].as3 * x[sl[i].iS3]; else sS3 = 0.0;
                if (sl[i].iT3 > -1) sT3 = sl[i].at3 * x[sl[i].iT3]; else sT3 = 0.0;
                if (sl[i].iB3 > -1) sB3 = sl[i].ab3 * x[sl[i].iB3]; else sB3 = 0.0;

                doublereal sE4, sW4, sN4, sS4, sT4, sB4;

                if (sl[i].iE4 > -1) sE4 = sl[i].ae4 * x[sl[i].iE4]; else sE4 = 0.0;
                if (sl[i].iW4 > -1) sW4 = sl[i].aw4 * x[sl[i].iW4]; else sW4 = 0.0;
                if (sl[i].iN4 > -1) sN4 = sl[i].an4 * x[sl[i].iN4]; else sN4 = 0.0;
                if (sl[i].iS4 > -1) sS4 = sl[i].as4 * x[sl[i].iS4]; else sS4 = 0.0;
                if (sl[i].iT4 > -1) sT4 = sl[i].at4 * x[sl[i].iT4]; else sT4 = 0.0;
                if (sl[i].iB4 > -1) sB4 = sl[i].ab4 * x[sl[i].iB4]; else sB4 = 0.0;

                const int it6 = sl[i].iP;

                long double residualQ = (sE + sW + sN + sS + sT + sB +
                    sE2 + sW2 + sN2 + sS2 + sT2 + sB2 +
                    sE3 + sW3 + sN3 + sS3 + sT3 + sB3 +
                    sE4 + sW4 + sN4 + sS4 + sT4 + sB4 +                               
                               b[it6]) - sl[i].ap * x[it6];

                

                const doublereal residual = static_cast<doublereal>(residualQ);

                // Чебышева
                dmax_loc += fabs((residual)); // сумма модулей
                dmax2_loc += residual * residual;
                x_new[it6] = x[it6] + rURF * (residual);


            }

            dmax += dmax_loc;
            dmax2 += dmax2_loc;

            dmax_loc = 0.0;
            dmax2_loc = 0.0;

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
            for (int i = 0; i < maxbound; ++i) {


                const int is = i * 2;
                const int it1 = col_ind_b[is];//iI
                const int it2 = col_ind_b[is + 1];//iW                   


                const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
                // Чебышева
                dmax_loc += fabs(residual); // сумма модулей.
                dmax2_loc += residual * residual;

                if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
                else x_new[it2] = x[it2] + rURF * (residual);


            }

            dmax += dmax_loc;
            dmax2 += dmax2_loc;

            //*std::min(1.0, 17000.0 / size0) 
            if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.

            if (bprintMessage) {
                if (iVar == PAM) {
                    //dmax/=maxelm;
                    if (j % degree == 0) {
#if doubleintprecision == 1
                        printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
                        printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

                    }
                }
            }

#pragma omp parallel for
            for (integer i = 0; i < size0; ++i)
            {
                x[i] = x_new[i];
            }

            if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
            delta = sqrt(dmax2 / (size0));

            if (bprintMessage) {
                if (j % degree == 0) {
                    std::cout << std::endl;
                }
            }

            j++;
        }// degree chebyshev


        delta /= r_start;



        // Корректировка нижней границы спектра lo линейного оператора А.
        // Жуков, Феодоритова. 13.02.2022
        correct_lo(lo, bmem, degree, hi, delta, bprintMessage);

       


        ++iout;
    }
   

    delete[] data_b;
    delete[] col_ind_b;

    std::cout << " " << lo << " " << j << " ";

    if (j > kend * degree - 100) {
        std::cout << "\n chebyshev divergence " << dmax << " ";
        system("pause");
    }

    delete[] x_new;

    if (bprintMessage) {
        if (iVar == PAM) {
            printf("calc complete...\n");
            system("pause");
        }
    }
    //printf("4000 %e \n", dmax);
    //system("pause");

} // chebyshev3Dnow_Alice

#endif


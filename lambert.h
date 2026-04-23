#pragma once

#include "doubleR.h"

inline std::pair<double, double> calc_c2c3(double psi, double eps = 1.e-6){
    double c2, c3, psi_sqrt = sqrt(fabs(psi));

    if (psi > eps) {
        c2 = (1. - cos(psi_sqrt)) / psi;
        c3 = (psi_sqrt - sin(psi_sqrt)) / (psi_sqrt * psi);
    }
    else if (psi < -eps){ 
        c2 = (1. - cosh(psi_sqrt)) / psi;
        c3 = (sinh(psi_sqrt) - psi_sqrt) / (psi_sqrt * abs(psi));
    }
    else {
        c2 = 0.5;
        c3 = 1. / 6.;
    }

    return {c2, c3};
}

inline std::tuple<double, double, double, double> Lambert ( 
                double mu,
                std::array <double, 3> r1_vec,
                std::array <double, 3> r2_vec,
                double dt,
                bool short_flag,
                double psi_min = -4. * M_PI,
                double psi_max = 4 * M_PI * M_PI,
                double dt_epsilon = 1.e-7) {

                    double t_m, r1, r2, cos_dnu, A, y,
                    psi = 0.0, mu_sqrt = sqrt(mu);
                    
                    if (short_flag == true) t_m = 1.0;
                    else t_m = -1.0;

                    r1 = sqrt(dot_prod(r1_vec, r1_vec));
                    r2 = sqrt(dot_prod(r2_vec, r2_vec));

                    cos_dnu = dot_prod(r1_vec, r2_vec) / (r1 * r2);

                    A = t_m * sqrt(r1 * r2 * (1. + cos_dnu));
                    if (A == 0.0) std::cout << "Немозможно решить задачу Ламберта." << std::endl;
                    auto [c2, c3] = calc_c2c3(psi);

                    while (true) {
                        y = r1 + r2 + A * (psi * c3 - 1.0) / sqrt(c2);

                        if (A > 0.0 && y < 0.0){
                            std::cout << "Предупреждение: y < 0" << std::endl;
                            psi_min += -y * sqrt(c2) / (A * c3);
                            continue;
                        }
                        double chi = sqrt(y / c2),
                        dt_f = (chi * chi * chi * c3 + A * sqrt(y)) / mu_sqrt;

                        if (fabs(dt - dt_f) < dt_epsilon) break;
                        
                        if (dt_f <= dt) psi_min = psi;
                        else psi_max = psi;

                        double psi_new = (psi_min + psi_max) * 0.5;
                        if (fabs(psi_new - psi) < 1.e-10) std::cout << "Не удалось найти решение" << std::endl;
                        psi = psi_new;
                        auto [c2, c3] = calc_c2c3(psi);
                    }
                    double f = 1. - y / r1, g = A * sqrt(y / mu), g_dot = 1. - y / r2;
                    std::array<double, 3> v1, v2;
                    for(size_t i = 0; i < 3; ++i){
                        v1[i] = (r2_vec[i] - f * r1_vec[i]) / g;
                        v2[i] = (g_dot * r2_vec[i] - r1_vec[i]) / g;
                    }

                    std::array<double, 6> state_v1, state_v2, elems1, elems2;
                    for (int i = 0; i < 3; ++i){
                        state_v1[i] = r1_vec[i];
                        state_v1[i+3] = v1[i];
                        state_v2[i] = r2_vec[i];
                        state_v2[i+3] = v2[i];
                    }
                    elems1 = get_elements(state_v1);
                    elems2 = get_elements(state_v2);

                    double p1 = elems1[0] * (1 - elems1[1] * elems1[1]),
                    p2 = elems2[0] * (1 - elems2[1] * elems2[1]),
                    v_t1 = sqrt(mu * p1) / r1,
                    v_t2 = sqrt(mu * p2) / r2,
                    E1 = solve_kepler(elems1[5], elems1[1]),
                    E2 = solve_kepler(elems2[5], elems2[1]),
                    v_r1 = sqrt(mu / p1) * elems1[1] * sqrt(1 - elems1[1] * elems1[1]) * sin(E1) / (1 - E1 * cos(E1)),
                    v_r2 = sqrt(mu / p2) * elems2[1] * sqrt(1 - elems2[1] * elems2[1]) * sin(E2) / (1 - E2 * cos(E2));

                    return {v_t1, v_r1, v_t2, v_r2};
}

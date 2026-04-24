#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <optional>
#include <iomanip>
#include <stdexcept>

constexpr double _RAD2GRAD = 180. / M_PI;
constexpr double mu = 398.6004415;


inline double dot_prod(const std::array<double, 3>& x, const std::array<double, 3>& y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline double solve_kepler(double M, double e){
    double E = M, dE = 1;
    while (dE > 1e-10){
        dE = (E - e * sin(E) - M) / (1 - e * cos(E));
        E -= dE;
    }
    return E;
}

inline std::array<double, 3> cross_prod(const std::array<double, 3>& x, const std::array<double, 3>& y){
    std::array<double, 3> cross;
    cross[0] = (x[1] * y[2] - x[2] * y[1]);
    cross[1] = (x[2] * y[0] - x[0] * y[2]);
    cross[2] = (x[0] * y[1] - x[1] * y[0]);

    return cross;
} 

inline std::array<double, 6> get_elements(std::array<double, 6> state_v){

    std::array <double, 3> r, v;
    std::copy((state_v).begin(), (state_v).begin() + 3, r.begin());
    std::copy((state_v).begin()+3, (state_v).end(), v.begin());

    double r_abs = sqrt(dot_prod(r, r));
    double ksi = dot_prod(v,v) / 2 - mu / r_abs;

    const std::array<double, 3> k = {0, 0, 1};

    std::array<double, 6> elements; //[a,e,i,OMEGA,omega,M]
    std::array<double, 3> h = cross_prod(r, v), n = cross_prod(k, h);
    std::array<double, 3> e_vec;

    double h_abs = sqrt(dot_prod(h,h)), n_abs = sqrt(dot_prod(n,n));

    e_vec[0] = ((dot_prod(v,v) - mu / r_abs) * r[0] - dot_prod(r, v) * v[0]) / mu;
    e_vec[1] = ((dot_prod(v,v) - mu / r_abs) * r[1] - dot_prod(r, v) * v[1]) / mu;
    e_vec[2] = ((dot_prod(v,v) - mu / r_abs) * r[2] - dot_prod(r, v) * v[2]) / mu;

    double e_abs = sqrt(dot_prod(e_vec, e_vec));
    elements[0] = -mu / 2 / ksi;
    elements[1] = e_abs;
    elements[2] = acos(h[2] / h_abs) * _RAD2GRAD; // i
    elements[3] = acos(n[0] / n_abs) * _RAD2GRAD; //OMEGA (ДВУ)
    if(n[1] < 0) elements[3] = 360 - elements[3];
    elements[4] = acos(dot_prod(n, e_vec) / n_abs / e_abs) * _RAD2GRAD; //omega (аргумент перицентра)
    if(e_vec[2] < 0) elements[4] = 360 - elements[4];

    double nu = acos(dot_prod(e_vec, r) / e_abs / r_abs);
    if (dot_prod(r, v) < 0) nu = 360 - nu;

    double E = 2 * atan(sqrt((1-e_abs) / (1+e_abs)) * tan(nu/2));
    elements[5] = (E - e_abs * sin(E)) * _RAD2GRAD;
    if (elements[5] < 0) elements[5] += 360;
    return elements;
}

inline std::array<double, 3> norm(const std::array<double, 3>& x, const std::array<double, 3>& y){
    std::array<double, 3> cross;
    double x_abs = sqrt(dot_prod(x, x)), y_abs = sqrt(dot_prod(y, y));
    cross[0] = (x[1] * y[2] - x[2] * y[1]) / x_abs / y_abs;
    cross[1] = (x[2] * y[0] - x[0] * y[2]) / x_abs / y_abs;
    cross[2] = (x[0] * y[1] - x[1] * y[0]) / x_abs / y_abs;
    double cross_abs = sqrt(dot_prod(cross, cross));
    for(size_t i = 0; i < 3; i++){cross[i] /= cross_abs;}

    return cross;
}

inline void print_vec(const std::array<double, 3>& x){
    std::cout << "[" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
}


struct Observation{
    double t;
    std::array<double, 3> R;
    std::array<double, 3> L;
};


struct OrbitSolution{
    double rho1, rho3;
    std::array<double, 3> r1, v1;
};


class DoubleRIteration{
    private:

        std::array<double, 3> L1, L2, L3, r_st1, r_st2, r_st3;
        double tau1, tau3;

        const double eps_dr = 1e-10;
        double p, a, e_sq;
        double cos_dE32, dE32, sin_dE32;

        double c1, c2;
        double rst1_abs, rst2_abs;
        double F1, F2;
        double dF1_dr1, dF2_dr1, dF1_dr2, dF2_dr2;

        std::array<double, 3> r2_vec, r3_vec, v2;
        std::array<double, 6> state_v;

    public:

        DoubleRIteration(std::array<double, 3> Ra, std::array<double, 3> Dec,
                         std::array<double, 3> r_st1, std::array<double, 3> r_st2,
                         std::array<double, 3> r_st3, double tau1, double tau3);

        DoubleRIteration& operator=(std::initializer_list<double> list);

        void calc_func(double r1, double r2);
        void calc_der(double r1, double r2);
        void solver(double r1, double r2);
        std::array<double, 6> get_state();
        void print_state();
        void print_elements();
};


class GoodingOD{

    private:

        std::array<double, 3> L1, L2, L3, r_st1, r_st2, r_st3;
        double tau1, tau3;

        const double eps_dr = 1e-12;
        double p, a, e_sq;
        double cos_dE32, dE32, sin_dE32;

        double c1, c2;
        double rst1_abs, rst2_abs;
        double F1, F2;
        double dF1_dr1, dF2_dr1, dF1_dr2, dF2_dr2;

        std::array<double, 3> r2_vec, r3_vec, v2;
        std::optional <std::array<double, 6>> state_v;

    public:

        GoodingOD(std::array<double, 3> Ra, std::array<double, 3> Dec, 
                            std::array<double, 3> r_st1, std::array<double, 3> r_st2, 
                            std::array<double, 3> r_st3, double tau1, double tau3);

        GoodingOD& operator=(std::initializer_list<double> list);

        void calc_func(double r1, double r2);
        void calc_der(double r1, double r2);
        void solver(double r1, double r2);
        std::array<double, 3> CALCPS(int nhrev, double rho1, double rho3);
        std::pair<double, double> OBS3LS(int nhrev, double rho1, double rho3);
        std::array<double, 6> get_state();
        void print_state();
        void print_elements();
};

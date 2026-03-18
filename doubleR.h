#pragma once
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <array>
#include <optional>

constexpr double mu = 398.6004415;
constexpr double _RAD2GRAD = 180 / M_PI;


class DoubleRIteration{
    private:

        std::array<double, 3> L1, L2, L3;
        std::array<double, 3> r_site1, r_site2, r_site3;
        double JD1, JD2, JD3;

        const double eps_dr = 1e-8;
        double p, a, e;
        double cos_dE32, dE32, sin_dE32;

        double tau1, tau3, c1, c2;
        double rst1_abs, rst2_abs;
        double F1, F2;
        double dF1_dr1, dF2_dr1, dF1_dr2, dF2_dr2;

        std::array<double, 3> r2_vec, r3_vec, v2;
        std::optional <std::array<double, 6>> state_v;

        double dot_prod(const std::array<double, 3>& x, const std::array<double, 3>& y) const {
            return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
        }

        std::array<double, 3> cross_prod(const std::array<double, 3>& x, const std::array<double, 3>& y)const{
            std::array<double, 3> cross;
            cross[0] = (x[1] * y[2] - x[2] * y[1]);
            cross[1] = (x[2] * y[0] - x[0] * y[2]);
            cross[2] = (x[0] * y[1] - x[1] * y[0]);

            return cross;
        } 

        std::array<double, 3> norm(const std::array<double, 3>& x, const std::array<double, 3>& y)const{
            std::array<double, 3> cross;
            double x_norm = sqrt(dot_prod(x, x)), y_norm = sqrt(dot_prod(y, y));
            cross[0] = (x[1] * y[2] - x[2] * y[1]) / x_norm / y_norm;
            cross[1] = (x[2] * y[0] - x[0] * y[2]) / x_norm / y_norm;
            cross[2] = (x[0] * y[1] - x[1] * y[0]) / x_norm / y_norm;
            double cross_norm = sqrt(dot_prod(cross, cross));
            for(size_t i = 0; i < 3; i++){cross[i] /= cross_norm;}
            return cross;
        }

    public:

        DoubleRIteration(std::array<double, 3> L1, 
                         std::array<double, 3> L2, 
                         std::array<double, 3> L3, 
                         double JD1, 
                         double JD2, 
                         double JD3, 
                         std::array<double, 3> r_site1, 
                         std::array<double, 3> r_site2, 
                         std::array<double, 3> r_site3
            ) {

                for(int i = 0; i < 3; i++){
                    this->L1[i] = L1[i];
                    this->L2[i] = L2[i];
                    this->L3[i] = L3[i];
                    this->r_site1[i] = r_site1[i];
                    this->r_site2[i] = r_site2[i];
                    this->r_site3[i] = r_site3[i];
                }

                this->JD1 = JD1;
                this->JD2 = JD2;
                this->JD3 = JD3;

                this->tau1 = JD1 - JD2;
                this->tau3 = JD3 - JD2;
                this->rst1_abs = sqrt(dot_prod(r_site1, r_site1));
                this->rst2_abs = sqrt(dot_prod(r_site2, r_site2));
                this->c1 = 2 * dot_prod(L1, r_site1);
                this->c2 = 2 * dot_prod(L2, r_site2);
                this->p = 0.0;
                this->a = 0.0;
                this->e = 0.0;
                this->dE32 = 0.0;
                this->cos_dE32 = 0.0;
                this->sin_dE32 = 0.0;
                this->F1 = 0.0;
                this->F2 = 0.0;
                this->dF1_dr1 = 0.0;
                this->dF1_dr2 = 0.0;
                this->dF2_dr1 = 0.0;
                this->dF2_dr2 = 0.0;
            }


        DoubleRIteration& operator=(std::initializer_list<double> list){
            if(list.size() == 6){
                std::array<double, 6> arr;
                std::copy(list.begin(), list.end(), arr.begin());
                state_v = arr;
            } else{
                throw std::invalid_argument("Требуется 6 элементов для вектора состояния.");
            }
            return *this;
        }

        void calc_func(double r1, double r2){

            double rho1, rho2, rho3, c3;
            double cos_dnu21, cos_dnu32;
            double sin_dnu21, sin_dnu31, sin_dnu32;
            double e_sin_nu2;
            double r1_abs, r2_abs, r3_abs;

            std::array<double, 3> W;
            std::array<double, 3> r1_vec;
            std::array<double, 3> n;

            try{
                double un_sq1 = c1 * c1 - 4 * (rst1_abs * rst1_abs - r1 * r1),
                un_sq2 = c2 * c2 - 4 * (rst2_abs * rst2_abs - r2 * r2);
                rho1 = (-c1 + sqrt(un_sq1)) / 2;
                rho2 = (-c2 + sqrt(un_sq2)) / 2;
                if (un_sq1 < 0 || un_sq2 < 0){
                    throw std::domain_error("Отрицательное значение под корнем для расстояния");
                    }
            }
            catch (const std::domain_error& error){
                std::cout << "c1 = " << c1 << ", "
                          << "rst1_abs = " << rst1_abs << ", "
                          << "r1 = " << r1 << std::endl;
                throw;
            }

            r1_vec[0] = rho1 * L1[0] + r_site1[0]; 
            r1_vec[1] = rho1 * L1[1] + r_site1[1];
            r1_vec[2] = rho1 * L1[2] + r_site1[2]; 
            r1_abs = sqrt(dot_prod(r1_vec, r1_vec));

            if (std::isnan(r1_abs)){std::exit(1);}

            r2_vec[0] = rho2 * L2[0] + r_site2[0]; 
            r2_vec[1] = rho2 * L2[1] + r_site2[1];
            r2_vec[2] = rho2 * L2[2] + r_site2[2];
            r2_abs = sqrt(dot_prod(r2_vec, r2_vec));

            W = cross_prod(r1_vec, r2_vec);
            for(size_t i = 0; i < W.size(); ++i) W[i] /= (r1_abs * r2_abs);
             
            rho3 = -dot_prod(r_site3, W) / dot_prod(L3, W);

            this->r3_vec[0] = rho3 * L3[0] + r_site3[0];
            this->r3_vec[1] = rho3 * L3[1] + r_site3[1];
            this->r3_vec[2] = rho3 * L3[2] + r_site3[2];
            r3_abs = sqrt(dot_prod(this->r3_vec, this->r3_vec));

            cos_dnu21 = dot_prod(r2_vec, r1_vec) / r2_abs / r1_abs;
            cos_dnu32 = dot_prod(this->r3_vec, this->r2_vec) / r3_abs / r2_abs;

            n = norm(r1_vec, r2_vec);

            sin_dnu21 = dot_prod(cross_prod(r1_vec, r2_vec), n) / r1_abs / r2_abs;
            sin_dnu31 = dot_prod(cross_prod(r1_vec, r3_vec), n) / r1_abs / r3_abs;
            sin_dnu32 = dot_prod(cross_prod(r2_vec, r3_vec), n) / r2_abs / r3_abs;

            if (sin_dnu31 < 0){
                c1 = r2_abs * sin_dnu32 / r1_abs / sin_dnu31;
                c3 = r2_abs * sin_dnu21 / r3_abs / sin_dnu31;
                p = (c1 * r1_abs + c3 * r3_abs - r2_abs) / (c1 + c3 - 1);
            }
            else {
                c1 = r1_abs * sin_dnu31 / r2_abs / sin_dnu32;
                c3 = r1_abs * sin_dnu21 / r3_abs / sin_dnu32;
                p = (c3 * r3_abs - c1 * r2_abs + r1_abs) / (-c1 + c3 + 1);
            }

            double e_cos_nu1 = p / r1_abs - 1;
            double e_cos_nu2 = p / r2_abs - 1;
            double e_cos_nu3 = p / r3_abs - 1;

            if (fabs(cos_dnu21 + 1) < 1e-6)
                e_sin_nu2 = (cos_dnu32 * e_cos_nu2 - e_cos_nu3) / sin_dnu31;
            else
                e_sin_nu2 = (-cos_dnu21 * e_cos_nu2 + e_cos_nu1) / sin_dnu21;

            e = sqrt(e_cos_nu2 * e_cos_nu2 + e_sin_nu2 * e_sin_nu2);
            a = p / (1 - e * e);
            double mean_anomaly = sqrt(mu / a / a / a);
            double S = r2_abs * sqrt(1 - e * e) * e_sin_nu2 / p;
            double C = r2_abs * (e * e + e_cos_nu2) / p;

            double sin_dE21 = r1_abs * (sin_dnu21 / sqrt(a * p) + (1 - cos_dnu21) * S / p);
            double cos_dE21 = 1 - r2_abs * r1_abs * (1 - cos_dnu21) / a / p;
            double dE21 = atan2(sin_dE21, cos_dE21);

            this->sin_dE32 = r3_abs * (sin_dnu32 / sqrt(a * p) - (1 - cos_dnu32) * S / p);
            this->cos_dE32 = 1 - r2_abs * r3_abs * (1 - cos_dnu32) / a / p;

            this->dE32 = atan2(sin_dE32, cos_dE32);
            double dM12 = -dE21 + 2 * S * sin(dE21 / 2) * sin(dE21 / 2) + C * sin_dE21;
            double dM32 = dE32 + 2 * S * sin(dE32 / 2) * sin(dE32 / 2) - C * sin_dE32;

            this->F1 = this->tau1 - dM12 / mean_anomaly;
            this->F2 = this->tau3 - dM32 / mean_anomaly;
        }

        void calc_der(double r1, double r2){
            
            double c1_cur = this->c1, eps = 1e-5;

            calc_func(r1 + eps * r1, r2);
            double F1r1_r = this->F1, F2r1_r = this->F2;

            this->c1 = c1_cur;
            calc_func(r1 - eps * r1, r2);
            double F1r1_l = this->F1, F2r1_l = this->F2;

            this->c1 = c1_cur;
            calc_func(r1, r2 + eps * r2);
            double F1r2_r = this->F1, F2r2_r = this->F2;

            this->c1 = c1_cur;
            calc_func(r1, r2 - eps * r2);
            double F1r2_l = this->F1, F2r2_l = this->F2;

            this->dF1_dr1 = (F1r1_r - F1r1_l) / 2 / eps / r1;
            this->dF2_dr1 = (F2r1_r - F2r1_l) / 2 / eps / r1;
            this->dF1_dr2 = (F1r2_r - F1r2_l) / 2 / eps / r2;
            this->dF2_dr2 = (F2r2_r - F2r2_l) / 2 / eps / r2;
        }

        void solver(double r1, double r2) {
            double dr1 = 1, dr2 = 1;
            double f,g;
            double F_1, F_2;

            while (fabs(dr1) > eps_dr || fabs(dr2) > eps_dr) { 
                double c1_cur;
                calc_func(r1, r2);
                F_1 = F1;
                F_2 = F2;

                c1_cur = this->c1;
                calc_der(r1, r2);

                double  D  = this->dF1_dr1 * this->dF2_dr2 - this->dF2_dr1 * this->dF1_dr2;
                double  d1 = this->dF2_dr2 * F_1 - this->dF1_dr2 * F_2, 
                        d2 = this->dF1_dr1 * F_2 - this->dF2_dr1 * F_1;

                dr1 = -d1 / D;
                dr2 = -d2 / D;
                if (fabs(dr1) > 0.05 * r1) dr1 = 0.1 * r1 * (dr1 > 0 ? 1 : -1);
                if (fabs(dr2) > 0.05 * r2) dr2 = 0.1 * r2 * (dr2 > 0 ? 1 : -1);
                r1 += dr1;
                r2 += dr2;
                std::cout << r2 << std::endl;
                this->c1 = c1_cur;
            }
            calc_func(r1, r2);

            f = 1 - a * (1 - cos_dE32) / r2;
            g = tau3 - sqrt(a * a * a / mu) * (dE32 - sin_dE32); 

            v2[0] = (r3_vec[0] - f * r2_vec[0]) / g; 
            v2[1] = (r3_vec[1] - f * r2_vec[1]) / g; 
            v2[2] = (r3_vec[2] - f * r2_vec[2]) / g; 

            for(int i = 0; i < 3; i++){
                (*state_v)[i] = r2_vec[i];
                (*state_v)[i+3] = v2[i];
            }
        }

        std::array<double, 6> get_state(){return *state_v;}

        std::array<double, 6> get_elements(){
            std::array <double, 3> r, v;
            std::copy((*state_v).begin(), (*state_v).begin() + 3, r.begin());
            std::copy((*state_v).begin()+3, (*state_v).end(), v.begin());

            double r_abs = sqrt(dot_prod(r, r));
            double ksi = dot_prod(v,v) / 2 - mu / r_abs;

            const std::array<double, 3> k = {0, 0, 1};
            std::array<double, 6> elements; //[a,e,i,omega,OMEGA,M]
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
            elements[4] = acos(n[0] / n_abs) * _RAD2GRAD; //OMEGA
            if(n[1] < 0) elements[4] = 360 - elements[4];
            elements[3] = acos(dot_prod(n, e_vec) / n_abs / e_abs) * _RAD2GRAD; //omega
            if(e_vec[2] < 0) elements[3] = 360 - elements[3];

            double nu = acos(dot_prod(e_vec, r) / e_abs / r_abs);
            if (dot_prod(r, v) < 0) nu = 360 - nu;

            double E = 2 * atan(sqrt((1-e_abs) / (1+e_abs)) * tan(nu/2));
            elements[5] = (E - e_abs * sin(E)) *  _RAD2GRAD;
            return elements;
        }
        
        void print_state() {
            const auto& state = get_state();
            std::cout << "[";
            for(size_t i = 0; i < state.size(); ++i){
                std::cout << state[i] << (i < state.size()-1 ? ", " : "");
            }
            std::cout << "]" << std::endl;
        }

        void print_elements(){
            const auto& elements = get_elements();
            std::cout << "[";
            for(size_t i = 0; i < elements.size(); ++i){
                std::cout << elements[i] << (i < elements.size()-1 ? ", " : "");
            }
            std::cout << "]" << std::endl;
        }

};

#include "doubleR.h"
#include "lambert.h"

DoubleRIteration::DoubleRIteration(std::array<double, 3> Ra, std::array<double, 3> Dec, 
                    std::array<double, 3> r_st1, std::array<double, 3> r_st2, 
                    std::array<double, 3> r_st3, double tau1, double tau3) {

        this->L1[0] = cos(Dec[0]) * cos(Ra[0]);
        this->L1[1] = cos(Dec[0]) * sin(Ra[0]);
        this->L1[2] = sin(Dec[0]);

        this->L2[0] = cos(Dec[1]) * cos(Ra[1]);
        this->L2[1] = cos(Dec[1]) * sin(Ra[1]);
        this->L2[2] = sin(Dec[1]);

        this->L3[0] = cos(Dec[2]) * cos(Ra[2]);
        this->L3[1] = cos(Dec[2]) * sin(Ra[2]);
        this->L3[2] = sin(Dec[2]);

        for(int i = 0; i < 3; i++){
            this->r_st1[i] = r_st1[i];
            this->r_st2[i] = r_st2[i];
            this->r_st3[i] = r_st3[i];
            this->r2_vec[i] = 0;
            this->r3_vec[i] = 0;
            this->v2[i] = 0;
            this->state_v[i] = 0;
            this->state_v[i+3] = 0;
        }

        this->tau1 = tau1;
        this->tau3 = tau3;
        this->rst1_abs = sqrt(dot_prod(r_st1, r_st1));
        this->rst2_abs = sqrt(dot_prod(r_st2, r_st2));
        this->c1 = 2 * dot_prod(L1, r_st1);
        this->c2 = 2 * dot_prod(L2, r_st2);
        this->p = 0.0;
        this->a = 0.0;
        this->e_sq = 0.0;
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


DoubleRIteration& DoubleRIteration::operator=(std::initializer_list<double> list){
    if(list.size() == 6){
        std::array<double, 6> arr;
        std::copy(list.begin(), list.end(), arr.begin());
        state_v = arr;
    } else{
        throw std::invalid_argument("Требуется 6 элементов для вектора состояния.");
    }
    return *this;
}

void DoubleRIteration::calc_func(double r1, double r2){

    double rho1, rho2, rho3, c3;
    double cos_dnu21, cos_dnu32, cos_dnu31;
    double sin_dnu21, sin_dnu31, sin_dnu32;
    double e_sin_nu2;
    double r1_abs, r2_abs, r3_abs;
    double mean_motion, S, C, dM32, dM12;

    std::array<double, 3> W;
    std::array<double, 3> r1_vec;
    std::array<double, 3> n;

    double un_sq1 = c1 * c1 - 4 * (rst1_abs * rst1_abs - r1 * r1);
    double un_sq2 = c2 * c2 - 4 * (rst2_abs * rst2_abs - r2 * r2);
    rho1 = (-c1 + sqrt(un_sq1)) / 2;
    rho2 = (-c2 + sqrt(un_sq2)) / 2;

    for(size_t i = 0; i < 3; ++i) {
        r1_vec[i] = rho1 * L1[i] + r_st1[i]; 
        r2_vec[i] = rho2 * L2[i] + r_st2[i];
    }
    r1_abs = sqrt(dot_prod(r1_vec, r1_vec));
    r2_abs = sqrt(dot_prod(r2_vec, r2_vec));

    W = cross_prod(r1_vec, r2_vec);
    for(size_t i = 0; i < W.size(); ++i) W[i] /= (r1 * r2);

    rho3 = -dot_prod(r_st3, W) / dot_prod(L3, W);

    if (rho1 < 0 || rho2 < 0 || rho3 < 0){
        std::cout << "rho1 = " << rho1
                  << " rho2 = " << rho2
                  << " rho3 = " << rho3 << std::endl;
        throw std::domain_error("Дальность отрицательная");
    }

    for(size_t i = 0; i < 3; ++i) r3_vec[i] = rho3 * L3[i] + r_st3[i];
    r3_abs = sqrt(dot_prod(r3_vec, r3_vec));

    cos_dnu21 = dot_prod(r2_vec, r1_vec) / r2_abs / r1_abs;
    cos_dnu31 = dot_prod(r3_vec, r1_vec) / r3_abs / r1_abs;
    cos_dnu32 = dot_prod(r3_vec, r2_vec) / r3_abs / r2_abs;

    sin_dnu21 = sqrt(1 - cos_dnu21 * cos_dnu21);
    sin_dnu31 = sqrt(1 - cos_dnu31 * cos_dnu31);
    sin_dnu32 = sqrt(1 - cos_dnu32 * cos_dnu32);

    //n = norm(r1_vec, r2_vec);

    //sin_dnu21 = dot_prod(cross_prod(r1_vec, r2_vec), n) / r1_abs / r2_abs;
    //sin_dnu31 = dot_prod(cross_prod(r1_vec, r3_vec), n) / r1_abs / r3_abs;
    //sin_dnu32 = dot_prod(cross_prod(r2_vec, r3_vec), n) / r2_abs / r3_abs;

    if (sin_dnu31 < 0){
        c1 = r2_abs * sin_dnu32 / r1_abs / sin_dnu31;
        c3 = r2_abs * sin_dnu21 / r3_abs / sin_dnu31;
        p = (c1 * r1_abs + c3 * r3_abs - r2_abs) / (c1 + c3 - 1.0);
    }
    else {
        c1 = r1_abs * sin_dnu31 / r2_abs / sin_dnu32;
        c3 = r1_abs * sin_dnu21 / r3_abs / sin_dnu32;
        p = (c3 * r3_abs - c1 * r2_abs + r1_abs) / (1.0 + c3 - c1);
    }

    double e_cos_nu1 = p / r1_abs - 1;
    double e_cos_nu2 = p / r2_abs - 1;
    double e_cos_nu3 = p / r3_abs - 1;

    if (fabs(cos_dnu21 + 1) < 1e-6)
        e_sin_nu2 = (cos_dnu32 * e_cos_nu2 - e_cos_nu3) / sin_dnu31;
    else
        e_sin_nu2 = (-cos_dnu21 * e_cos_nu2 + e_cos_nu1) / sin_dnu21;

    e_sq = e_cos_nu2 * e_cos_nu2 + e_sin_nu2 * e_sin_nu2;
    a = p / (1 - e_sq);

    if (e_sq - 1 > 1e-10){
        mean_motion = sqrt(- mu / a / a / a);
        S = r2_abs * sqrt(e_sq - 1) * e_sin_nu2 / p;
        C = r2_abs * (e_sq + e_cos_nu2) / p;
        double sh_dF32 = r3_abs * (sin_dnu32 / sqrt(-a * p) - (1 - cos_dnu32) * S / p);
        double ch_dF32 = 1 - r3_abs * r2_abs * (1 - cos_dnu32) / a / p;
        double sh_dF21 = r1_abs * (sin_dnu21 / sqrt(-a * p) + (1 - cos_dnu21) * S / p);
        double ch_dF21 = 1 - r2_abs * r1_abs * (1 - cos_dnu21) / a / p;

        double dF32 = log(sh_dF32 + sqrt(sh_dF32 * sh_dF32 + 1));
        double dF21 = log(sh_dF21 + sqrt(sh_dF21 * sh_dF21 + 1));

        dM32 = -dF32 + S * (ch_dF32 - 1) + C * sh_dF32;
        dM12 =  dF21 + S * (ch_dF21 - 1) - C * sh_dF21;
    }
    else{
        mean_motion = sqrt(mu / a / a / a);
        S = r2_abs * sqrt(1 - e_sq) * e_sin_nu2 / p;
        C = r2_abs * (e_sq + e_cos_nu2) / p;

        double sin_dE21 = r1_abs * (sin_dnu21 / sqrt(a * p) + (1 - cos_dnu21) * S / p);
        double cos_dE21 = 1 - r2_abs * r1_abs * (1 - cos_dnu21) / a / p;
        double dE21 = atan2(sin_dE21, cos_dE21);
        sin_dE32 = r3_abs * (sin_dnu32 / sqrt(a * p) - (1 - cos_dnu32) * S / p);
        cos_dE32 = 1 - r2_abs * r3_abs * (1 - cos_dnu32) / a / p;
        dE32 = atan2(sin_dE32, cos_dE32);

        dM12 = -dE21 + S * (1 - cos_dE21) + C * sin_dE21;
        dM32 =  dE32 + S * (1 - cos_dE32) - C * sin_dE32;
    }
    F1 = tau1 - dM12 / mean_motion;
    F2 = tau3 - dM32 / mean_motion;
}

void DoubleRIteration::calc_der(double r1, double r2){
    
    double c1_cur = c1, eps = 1e-5;
    calc_func(r1 + eps * r1, r2);
    double F1r1_r = F1, F2r1_r = F2;

    c1 = c1_cur;
    calc_func(r1 - eps * r1, r2);
    double F1r1_l = F1, F2r1_l = F2;

    c1 = c1_cur;
    calc_func(r1, r2 + eps * r2);
    double F1r2_r = F1, F2r2_r = F2;

    c1 = c1_cur;
    calc_func(r1, r2 - eps * r2);
    double F1r2_l = F1, F2r2_l = F2;

    dF1_dr1 = (F1r1_r - F1r1_l) / 2 / eps / r1;
    dF2_dr1 = (F2r1_r - F2r1_l) / 2 / eps / r1;
    dF1_dr2 = (F1r2_r - F1r2_l) / 2 / eps / r2;
    dF2_dr2 = (F2r2_r - F2r2_l) / 2 / eps / r2;
}

void DoubleRIteration::solver(double r1, double r2) {
    double dr1 = 1, dr2 = 1;
    double f,g;
    double F_1, F_2;
    double k = 0.2;

    while (fabs(dr1) > eps_dr || fabs(dr2) > eps_dr) { 
        double c1_cur = c1;
        try{calc_func(r1, r2);}
        catch(const std::exception& e){
            std::cout << e.what() << std::endl; 
            return;
        }
        F_1 = F1;
        F_2 = F2;
        
        c1 = c1_cur;
        calc_der(r1, r2);

        double  D  = dF1_dr1 * dF2_dr2 - dF2_dr1 * dF1_dr2;
        double  d1 = dF2_dr2 * F_1 - dF1_dr2 * F_2, 
                d2 = dF1_dr1 * F_2 - dF2_dr1 * F_1;
        dr1 = -d1 / D;
        dr2 = -d2 / D;

        r1 += k * dr1;
        r2 += k * dr2;
        c1 = c1_cur;
    }
    calc_func(r1, r2);
    
    double r2_abs = sqrt(dot_prod(r2_vec, r2_vec));

    f = 1 - a * (1 - cos_dE32) / r2_abs;
    g = tau3 - sqrt(a * a * a / mu) * (dE32 - sin_dE32); 
    //std::cout << "a = " << a << " cos_dE32 = " << cos_dE32 << " r2_abs = " << r2_abs 
    //          << " dE32 = " << dE32 << " sin_dE32 = " << sin_dE32 << std::endl;
    for(size_t i = 0; i < 3; ++i) v2[i] = (r3_vec[i] - f * r2_vec[i]) / g;
    for(size_t i = 0; i < 3; i++){
        state_v[i] = r2_vec[i];
        state_v[i+3] = v2[i];
    }
}

std::array<double, 6> DoubleRIteration::get_state(){return state_v;}

void DoubleRIteration::print_state() {
    const auto& state = get_state();
    std::cout << "[";
    for(size_t i = 0; i < state.size(); ++i){
        std::cout << std::fixed << std::setprecision(14) << state[i] << (i < state.size()-1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
}

void DoubleRIteration::print_elements(){
    const auto& elements = get_elements(state_v);
    std::cout << "[";
    for(size_t i = 0; i < elements.size(); ++i){
        std::cout << std::fixed << std::setprecision(14) << elements[i] << (i < elements.size()-1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
}


GoodingOD::GoodingOD(std::array<double, 3> Ra, std::array<double, 3> Dec,
                    std::array<double, 3> r_st1, std::array<double, 3> r_st2, 
                    std::array<double, 3> r_st3, double tau1, double tau3) {

        this->L1[0] = cos(Dec[0]) * cos(Ra[0]);
        this->L1[1] = cos(Dec[0]) * sin(Ra[0]);
        this->L1[2] = sin(Dec[0]);

        this->L2[0] = cos(Dec[1]) * cos(Ra[1]);
        this->L2[1] = cos(Dec[1]) * sin(Ra[1]);
        this->L2[2] = sin(Dec[1]);

        this->L3[0] = cos(Dec[2]) * cos(Ra[2]);
        this->L3[1] = cos(Dec[2]) * sin(Ra[2]);
        this->L3[2] = sin(Dec[2]);

        for(int i = 0; i < 3; i++){
            this->r_st1[i] = r_st1[i];
            this->r_st2[i] = r_st2[i];
            this->r_st3[i] = r_st3[i];
            this->r2_vec[i] = 0;
            this->r3_vec[i] = 0;
            this->v2[i] = 0;
        }

        this->tau1 = tau1;
        this->tau3 = tau3;
        this->rst1_abs = sqrt(dot_prod(r_st1, r_st1));
        this->rst2_abs = sqrt(dot_prod(r_st2, r_st2));
        this->c1 = 2 * dot_prod(L1, r_st1);
        this->c2 = 2 * dot_prod(L2, r_st2);
        this->p = 0.0;
        this->a = 0.0;
        this->e_sq = 0.0;
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


GoodingOD& GoodingOD::operator=(std::initializer_list<double> list){
    if(list.size() == 6){
        std::array<double, 6> arr;
        std::copy(list.begin(), list.end(), arr.begin());
        state_v = arr;
    } else{
        throw std::invalid_argument("Требуется 6 элементов для вектора состояния.");
    }
    return *this;
}


void GoodingOD::solver(const Observation& obs1, const Observation& obs2, const Observation& obs3,
                        double rho1_init, double rho3_init, Vector& r1_out, Vector& v1_out){
    double rho1 = rho1_init, rho3 = rho3_init;
    
    for(int k = 0; k < 50; k++){
        Vector F;
        double conv = sqrt(F.x * F.x + F.y * F.y);

        if (conv < 1e-10){
            
            std::array<double, 3> r1, r3;

            for(size_t i; i < 3; ++i){
                r1[i] = obs1.R[i] + obs1.L[i] * rho1;
                r3[i] = obs3.R[i] + obs3.L[i] * rho3;
            }
            Vector v1, v3;

            double tau3 = obs3.t - obs1.t;
            auto [r1_out, v1_out] = Lambert(r1, r3, tau3, true);
        }
    }
    
}

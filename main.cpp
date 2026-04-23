#include "doubleR.h"

int main(){
    std::array<double, 3> Ra = {0.939913 / _RAD2GRAD, 
                                45.025748 / _RAD2GRAD,
                                67.886655 / _RAD2GRAD};
    std::array<double, 3> Dec = {18.667717 / _RAD2GRAD, 
                                35.664741 / _RAD2GRAD,
                                36.996583  / _RAD2GRAD};

    double tau1 = -480 * 1e-3, tau3 = 240 * 1e-3;
    std::array<double, 3> r_site1 = {4.054881, 2.748195, 4.074237};
    std::array<double, 3> r_site2 = {3.965224, 2.888232, 4.074364};
    std::array<double, 3> r_site3 = {3.905073, 2.956935, 4.07443};
    double r1 = 10, r2 = 11.020;

    DoubleRIteration orbit(Ra, Dec, r_site1, r_site2, r_site3, tau1, tau3);
    orbit.solver(r1, r2);
    orbit.print_state();
    orbit.print_elements();
    return 0;
}
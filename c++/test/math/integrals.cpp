
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Test for physim/c++/include/math
// last updated:    03/07/2022


#include "../../include/math/integral.h"


int main() {


    Integral integral; 

    Cosine cos(5., 2.); 
    Exponential exp(M_E);
    Sine sin(2., 1.); 
    Logarithm log(M_E, 1., 1.);

    Functor g('/', sin, log); 
    Functor f('c', exp, g);
    Functor h('*', cos, f); 
    
    std::cout << "\n\nMidpoint" << std::endl; 
    integral.midpoint_fixed(3., 7., h, 1.e-7); 
    integral.print_integral(1.e-7);
    integral.print_error(); 
    
    std::cout << "\n\nTrapexoid" << std::endl; 
    integral.trapexoid_fixed(3., 7., h, 1.e-7); 
    integral.print_integral(1.e-7);
    integral.print_error(); 

    std::cout << "\n\nSimpson" << std::endl; 
    integral.simpson_fixed(3., 7., h, 1.e-7); 
    integral.print_integral(1.e-7);
    integral.print_error(); 
    
    // std::cout << "\n\nMean" << std::endl; 
    // integral.mean_fixed(3., 7., h, 1.e-7); 
    // integral.print_integral(1.e-7);
    // integral.print_error(); 
    
    
    return 0; 

}
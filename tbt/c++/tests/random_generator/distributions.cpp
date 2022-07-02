
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Particle(class) is an object class that stores information about particles.
// last updated:    17/06/2022

#include "../../include/math/random_generator.h"
#include "../../include/gplot++.h"


int main(int argc, char** argv) {    

    unsigned int nmax = atoi(argv[1]); 
    RandomGenerator rg(1);
    std::vector<double> rand{}, exp{}, gaussBM{}, gaussAR{}; 
    Gnuplot plot{};
    plot.redirect_to_png("images/distributions.png");
    plot.multiplot(2, 2, "Distributions");


    // Rand histo
    for (unsigned int i{}; i < nmax; i++) 
        rand.push_back(rg.rand());     
    
    plot.set_xlabel("Value");
    plot.set_ylabel("Number of counts");
    plot.set_yrange(0., NAN); 
    plot.histogram(rand, 20, "Rand");
    plot.show();

    std::cout << "Plot ended (1/4)" << std::endl; 

    // Exp histo
    rg.set_seed(1); 
    for (unsigned int i{}; i < nmax; i++) 
        exp.push_back(rg.exp(1));     
    
    plot.set_xlabel("Value");
    plot.set_ylabel("Number of counts");
    plot.set_xrange(0., 10.); 
    plot.histogram(exp, 100, "Exp");
    plot.show();

    std::cout << "Plot ended (2/4)" << std::endl; 

    // GaussBM histo
    rg.set_seed(1); 
    for (unsigned int i{}; i < nmax; i++) 
        gaussBM.push_back(rg.gauss_box_muller(2., 1.));     
    
    plot.set_xlabel("Value");
    plot.set_ylabel("Number of counts");
    plot.histogram(gaussBM, 40, "Gauss BM");
    plot.show();

    std::cout << "Plot ended (3/4)" << std::endl; 

    // // GaussAR histo
    // rg.set_seed(1); 
    // for (unsigned int i{}; i < nmax; i++) 
    //     gaussAR.push_back(rg.gauss_accept_reject(2., 1.)); 
    
    // plot.set_xlabel("Value");
    // plot.set_ylabel("Number of counts");
    // plot.histogram(gaussAR, 40, "Gauss AR");
    // plot.show();

    // std::cout << "Plot ended (4/4)" << std::endl; 

    plot.show();
    
    return 0; 
}
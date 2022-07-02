
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of an harmonic oscillator using ODE's methods.
// last updated:    02/07/2022


#include "../../include/physics/oscillators.h"
#include "../gplot++.h"
#include <chrono>


int main() {

    Gnuplot plot{};
    plot.redirect_to_png("images/harmonic.png");
    plot.set_xlabel("Time [s]");
    plot.set_ylabel("Position [m]"); 
    
    HarmonicOscillator oscill(1); 
    Position pos({0., 0., 0.}, {1., 0., 0.});
       
    const double tmax{70}, h{0.001};
    std::vector<double> coord_x{}, time{};
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::time_t end_time;
    
    // Euler's method
    start = std::chrono::system_clock::now();
    while (oscill.get_time() < tmax) {
        time.push_back(oscill.get_time()); 
        coord_x.push_back(pos.get_coord_x());
        pos.set_position(oscill.euler(pos.get_position(), h)); 
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t(end);
    
    std::cout << "Simulation with Euler's method ended\n"; 
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n" << std::endl;
    
    plot.plot(time, coord_x, "Euler"); 
    pos.set_position({0., 0., 0.}, {1., 0., 0.}); 
    oscill.reset_time(); 
    time.clear(); 
    coord_x.clear(); 
    
    // Runge Kutta's method
    start = std::chrono::system_clock::now();
    while (oscill.get_time() < tmax) {
        time.push_back(oscill.get_time()); 
        coord_x.push_back(pos.get_coord_x());
        pos.set_position(oscill.euler(pos.get_position(), h)); 
    }    
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t(end);
    
    std::cout << "Simulation with Runge Kutta's method ended\n"; 
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n" << std::endl;
    
    plot.plot(time, coord_x, "RK4"); 
    plot.show(); 
    std::cout << "Plot ended" << std::endl; 

    return 0; 

}

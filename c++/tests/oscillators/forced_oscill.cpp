
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of a forced oscillator using ODE's methods.
// last updated:    19/06/2022

#include "../include/physics/oscillators.h"
#include "../include/gplot++.h"


int main() {

    Gnuplot plot{};
    plot.redirect_to_png("images/oscillators/forced.png");
    plot.set_xlabel("Time [s]");
    plot.set_ylabel("Position [m]"); 

    ForcedOscillator oscill({0., 0., 0.}, {0., 0., 0.}, 10., 5.); 

    const double tmax{100}, h{0.001};  
    std::vector<double> coord_x{}, time{};
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::time_t end_time;

    // Runge Kutta's method
    start = std::chrono::system_clock::now();
    while (oscill.get_time() < tmax) {
        time.push_back(oscill.get_time()); 
        coord_x.push_back(oscill.get_coord_x());
        oscill.rk4(h); 
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Simulation with Runge Kutta's method ended\n"; 
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n" << std::endl;

    plot.plot(time, coord_x, "RK4"); 
    plot.show(); 
    std::cout << "Plot ended\n"; 

    return 0; 

} 


// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of Earth and Mars' orbirts around the sun using ODE's methods.
// last updated:    05/07/2022


#include "../../include/physics/gravitational_field.h"
#include "../../include/physics/planet.h"
#include "../gplot++.h"
#include <chrono>


int main() {

    Planet sun("Sun"), earth("Earth"), mars("Mars");
    sun.set_position({0., 0., 0.}, {0., 0., 0.});
    sun.print_body(); 
    sun.print_position();
    earth.set_position({-earth.get_coord_perihelion(), 0., 0.}, {0., -earth.get_vel_perihelion(), 0.});  
    earth.print_body(); 
    earth.print_position(); 
    mars.set_position({mars.get_coord_aphelion(), 0., 0.}, {0., mars.get_vel_aphelion(), 0.});  
    mars.print_body(); 
    mars.print_position(); 
    GravitationalField sun_gravity(sun.get_mass(), sun.get_coordinates()); 

    // number of seconds in 1 year: 31'536'000 s
    // number of seconds in 1 month: 2'628'000 s 
    // number of seconds in 1 week: 604'800 s
    // number of seconds in 1 day: 86'400 s  
    
    const double tmax{3.1554e7}, h{1};
    unsigned int count{};
    std::vector<double> coord_x_earth{}, coord_y_earth{}, coord_x_mars{}, coord_y_mars{};
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::time_t end_time;

    // Runge Kutta's method
    start = std::chrono::system_clock::now();
    while (sun_gravity.get_time() < tmax) {
        if (count % 10000 == 0) {
            coord_x_earth.push_back(earth.get_coord_x());
            coord_y_earth.push_back(earth.get_coord_y());            
            coord_x_mars.push_back(mars.get_coord_x());
            coord_y_mars.push_back(mars.get_coord_y());
            std::cout << count << std::endl;
        }
        earth.set_position(sun_gravity.rk4(earth.get_position(), h)); 
        mars.set_position(sun_gravity.rk4(mars.get_position(), h));   
        sun_gravity.increase_time(h);       
        count++;
    }    
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t(end);
    
    std::cout << "Simulation with Runge Kutta's method ended\n"; 
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n" << std::endl;
    std::cout << "Final position: " << std::endl;
    earth.print_position(); 

    Gnuplot plot{};
    plot.redirect_to_png("images/earth_mars_sun.png");
    plot.set_xlabel("x [km]");
    plot.set_ylabel("y [km]");     
    plot.plot(coord_x_earth, coord_y_earth, "Earth's orbit"); 
    plot.plot(coord_x_mars, coord_y_mars, "Mars' orbit"); 
    plot.show(); 
    std::cout << "Plot ended" << std::endl; 

    return 0; 

}

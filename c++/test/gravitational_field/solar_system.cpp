
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of Solar System planet's orbirts around the sun using ODE's methods.
// last updated:    05/07/2022


#include "../../include/physics/gravitational_field.h"
#include "../../include/physics/system.h"
#include "../gplot++.h"
#include <chrono>


int main() {

    std::cout << "Simulation of the orbits of the Solar System planets with Runge Kutta's method\n"; 

    Planet sun("Sun"), mercury("Mercury"), venus("Venus"), earth("Earth"), mars("Mars"), jupiter("Jupiter"), saturn("Saturn"), uranus("Uranus"), neptune("Neptune");
    sun.set_position({0., 0., 0.}, {0., 0., 0.});
    mercury.set_position({mercury.get_coord_perihelion(), 0., 0.}, {0., mercury.get_vel_perihelion(), 0.});  
    venus.set_position({venus.get_coord_perihelion(), 0., 0.}, {0., venus.get_vel_perihelion(), 0.});  
    earth.set_position({earth.get_coord_perihelion(), 0., 0.}, {0., earth.get_vel_perihelion(), 0.});  
    mars.set_position({mars.get_coord_perihelion(), 0., 0.}, {0., mars.get_vel_perihelion(), 0.});  
    jupiter.set_position({jupiter.get_coord_perihelion(), 0., 0.}, {0., jupiter.get_vel_perihelion(), 0.});  
    saturn.set_position({saturn.get_coord_perihelion(), 0., 0.}, {0., saturn.get_vel_perihelion(), 0.});  
    uranus.set_position({uranus.get_coord_perihelion(), 0., 0.}, {0., uranus.get_vel_perihelion(), 0.});  
    neptune.set_position({neptune.get_coord_perihelion(), 0., 0.}, {0., neptune.get_vel_perihelion(), 0.});  

    GravitationalField sun_gravity(sun); 
    System<Planet> solar_system; 
    solar_system.add_object(mercury); 
    solar_system.add_object(venus); 
    solar_system.add_object(earth); 
    solar_system.add_object(mars); 
    solar_system.add_object(jupiter); 
    solar_system.add_object(saturn); 
    solar_system.add_object(uranus); 
    solar_system.add_object(neptune); 
    std::cout << "Number of planets = " << solar_system.get_objects_count() << std::endl; 

    unsigned int days_to_seconds{86400}, count{};
    const double h{1.}; 
    std::vector<double> coord_x{sun.get_coord_x()}, coord_y{sun.get_coord_y()};
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::time_t end_time;

    Gnuplot plot{};
    plot.redirect_to_png("images/solar_system.png");
    plot.set_xlabel("x [km]");
    plot.set_ylabel("y [km]");  
    plot.plot(coord_x, coord_y, "Sun", Gnuplot::LineStyle::POINTS); 
    coord_x.clear(); 
    coord_y.clear();  

    for (auto i : solar_system.get_objects()) { 

        i.print_body(); 
        i.print_position(); 

        // Runge Kutta's method
        std::cout << "\nSimulation started" << std::endl; 
        start = std::chrono::system_clock::now();

        while (sun_gravity.get_time() < i.get_period() * days_to_seconds) {
            if (count % 10000 == 0) {
                coord_x.push_back(i.get_coord_x());
                coord_y.push_back(i.get_coord_y());            
            }
            i.set_position(sun_gravity.rk4(i.get_position(), h)); 
            sun_gravity.increase_time(h);       
            count++;
        } 

        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        end_time = std::chrono::system_clock::to_time_t(end);
        
        std::cout << "Simulation ended" << std::endl; 
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n" << std::endl;
        std::cout << "Final position: " << std::endl;
        i.print_position(); 

        std::cout << "Plot started" << std::endl; 
        plot.plot(coord_x, coord_y, i.get_name()); 
        plot.show(); 
        std::cout << "Plot ended" << std::endl; 

        sun_gravity.reset_time();
        count = 0; 
        coord_x.clear(); 
        coord_y.clear();
          
    }
    
    plot.show(); 
    return 0; 

}

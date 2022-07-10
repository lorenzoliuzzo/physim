
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of Solar System planet's orbirts around the sun using ODE's methods.
// last updated:    06/07/2022


#include "../../include/physics/tools/system.h"
#include "../../include/physics/objects/celestial_body.h"
#include "../gplot++.h"


int main() {

    std::cout << "Simulation of the orbits of the Solar System planets with Runge Kutta's method\n"; 

    Planet sun("Sun"), mercury("Mercury"), venus("Venus"), earth("Earth"), mars("Mars"), jupiter("Jupiter"), saturn("Saturn"), uranus("Uranus"), neptune("Neptune");
    mercury.set_position({mercury.get_coord_perihelion(), 0., 0.}, {0., mercury.get_vel_perihelion(), 0.});  
    venus.set_position({venus.get_coord_perihelion(), 0., 0.}, {0., venus.get_vel_perihelion(), 0.});  
    earth.set_position({earth.get_coord_perihelion(), 0., 0.}, {0., earth.get_vel_perihelion(), 0.});  
    mars.set_position({mars.get_coord_perihelion(), 0., 0.}, {0., mars.get_vel_perihelion(), 0.});  
    jupiter.set_position({jupiter.get_coord_perihelion(), 0., 0.}, {0., jupiter.get_vel_perihelion(), 0.});  
    saturn.set_position({saturn.get_coord_perihelion(), 0., 0.}, {0., saturn.get_vel_perihelion(), 0.});  
    uranus.set_position({uranus.get_coord_perihelion(), 0., 0.}, {0., uranus.get_vel_perihelion(), 0.});  
    neptune.set_position({neptune.get_coord_perihelion(), 0., 0.}, {0., neptune.get_vel_perihelion(), 0.});  

    System<Planet> solar_system; 
    solar_system.add_object(mercury); 
    solar_system.add_object(venus); 
    solar_system.add_object(earth); 
    solar_system.add_object(mars); 
    solar_system.add_object(jupiter); 
    solar_system.add_object(saturn); 
    solar_system.add_object(uranus); 
    solar_system.add_object(neptune); 
        
    for (auto i : solar_system.get_objects()) {
        i.print_body(); 
        i.activate_gravitational_field(); 
        for (auto j : solar_system.get_objects()) {
            i.add_gravitational_attraction(i.gravitational_attraction(i.Position::get_coordinates())); 
        }
        i.print_gravitational_attraction();
    }



    // unsigned int days_to_seconds{86400}, count{};
    // double h{10.}; 
    // std::vector<double> coord_x{sun.get_coord_x()}, coord_y{sun.get_coord_y()};
    // std::chrono::duration<double> elapsed_seconds;
    // std::chrono::time_point<std::chrono::system_clock> start, end;
    // std::time_t end_time;

    // Gnuplot plot{};
    // plot.redirect_to_png("images/solar_system.png");
    // plot.set_xlabel("x [km]");
    // plot.set_ylabel("y [km]");  
    // plot.plot(coord_x, coord_y, "Sun", Gnuplot::LineStyle::POINTS); 
    // coord_x.clear(); 
    // coord_y.clear();  
    
    // std::cout << "\nSimulation started" << std::endl; 
    // std::cout << "Number of planets = " << solar_system.get_objects_count() << std::endl; 
    // start = std::chrono::system_clock::now();

    // for (auto i : solar_system.get_objects()) { 

    //     i.print_body(); 
    //     i.print_position(); 

    //     // Runge Kutta's method
    //     while (sun_gravity.get_time() < i.get_period() * days_to_seconds) {
    //         if (count % 1000 == 0) {
    //             coord_x.push_back(i.get_coord_x());
    //             coord_y.push_back(i.get_coord_y());            
    //         }
    //         i.set_position(sun_gravity.rk4(i.get_position(), h)); 
    //         sun_gravity.increase_time(h);       
    //         count++;
    //     } 

    //     end = std::chrono::system_clock::now();
    //     elapsed_seconds = end - start;
    //     end_time = std::chrono::system_clock::to_time_t(end);
        
    //     std::cout << "\nSimulation ended" << std::endl; 
    //     std::cout << "elapsed time: " << elapsed_seconds.count() << " s" << std::endl;
    //     std::cout << "Final position: " << std::endl;
    //     i.print_position(); 

    //     std::cout << "Plot started" << std::endl; 
    //     plot.plot(coord_x, coord_y, i.get_name()); 
    //     std::cout << "Plot ended\n" << std::endl; 

    //     sun_gravity.reset_time();
    //     count = 0; 
    //     coord_x.clear(); 
    //     coord_y.clear();
    //     h = h * 2; 

    // }
    // end = std::chrono::system_clock::now();
    // elapsed_seconds = end - start;
    // end_time = std::chrono::system_clock::to_time_t(end);
    
    // std::cout << "\n\nSimulation ended" << std::endl; 
    // std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n" << std::endl;
    // plot.show(); 
    return 0; 

}

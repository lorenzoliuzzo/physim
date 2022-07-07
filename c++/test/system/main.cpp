


#include "../../include/physics/system.h"
#include "../../include/physics/gravitational_field.h"
#include "../../include/physics/celestial_body.h"
#include "../gplot++.h"


int main() {
    
    System<Planet> solar_system; 
    Planet sun("Sun"), earth("Earth"), mars("Mars");
    sun.set_position({0., 0., 0.}, {0., 0., 0.}); 
    earth.set_position({earth.get_coord_aphelion(), 0., 0.}, {0., earth.get_vel_aphelion(), 0.}); 
    mars.set_position({-mars.get_coord_perihelion(), 0., 0.}, {0., - mars.get_vel_perihelion(), 0.}); 

    solar_system.add_object(sun); 
    solar_system.add_object(earth); 
    solar_system.add_object(mars); 

    std::cout << "\nNumber of planets = " << solar_system.get_objects_count() << std::endl; 

    for (auto i : solar_system.get_objects()) {
        i.print_body(); 
        i.print_position(); 
    }

    std::vector<std::vector<double>> appo = zeros(3, 2); 
    const double h{5.}, tmax{1.e6}; 
    while (solar_system.get_time() < tmax) {
        unsigned int i{}; 
        for (auto k : solar_system.get_objects()) {
            GravitationalField gravity(k); 
            for (unsigned int j{}; j < solar_system.get_objects_count(); j++) {
                if (j != i) appo += gravity.rk4(solar_system.get_object(j).get_position(), h); 
            }

            k.set_position(appo); 
            solar_system.increase_time(h); 
            i++; 
            appo.clear(); 
        }
    }

    std::cout << "\nFinal result " << std::endl; 
    for (auto i : solar_system.get_objects()) {
        i.print_body(); 
        i.print_position(); 
    }    
    
    // solar_system.activate_gravitational_field(); 
    // while(solar_system.get_time() < tmax) {
    //     solar_system.evolve(); 
    // }

    // std::cout << "time = " << solar_system.get_time() << std::endl; 
    // for (auto i : solar_system.get_objects()) {
    //     i.print_body(); 
    //     i.print_position(); 
    // }

    return 0; 

}

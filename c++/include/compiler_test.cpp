// compile test

// #include "math/integral.h"
// #include "math/statistics.h"
// #include "physics/position.h"
#include "physics/system.h"
// #include "physics/planet.h"
// #include "physics/gravitational_field.h"


int main() {

    const double tmax{10}; 
    
    System<Planet> SolarSystem; 
    Planet sun("Sun"), earth("Earth"), mars("Mars");
    sun.set_position({0., 0., 0.}, {0., 0., 0.}); 
    earth.set_position({earth.get_coord_aphelion(), 0., 0.}, {0., earth.get_vel_aphelion(), 0.}); 
    mars.set_position({-mars.get_coord_perihelion(), 0., 0.}, {0., - mars.get_vel_perihelion(), 0.}); 

    SolarSystem.add_object(sun); 
    SolarSystem.add_object(earth); 
    SolarSystem.add_object(mars); 

    for (auto i : SolarSystem.get_objects()) {
        i.print_body(); 
        i.print_position(); 
    }
    
    std::cout << "Number of planets = " << SolarSystem.get_objects_count() << std::endl; 
    SolarSystem.activate_gravitational_field(); 
    while(SolarSystem.get_time() < tmax) {
        SolarSystem.evolve(); 
    }

    std::cout << "time = " << SolarSystem.get_time() << std::endl; 
    for (auto i : SolarSystem.get_objects()) {
        i.print_body(); 
        i.print_position(); 
    }

    return 0; 
}
// compile test

// #include "math/integral.h"
// #include "math/statistics.h"
// #include "physics/position.h"
#include "physics/system.h"
#include "physics/planets.h"
#include "physics/gravitational_field.h"


int main() {
    
    System<Planet> SolarSystem; 
    Planet sun("Sun"), earth("Earth"), mars("Mars");
    sun.set_position({0., 0., 0.}, {0., 0., 0.}); 
    earth.set_position({earth.get_coord_afelio(), 0., 0.}, {earth.get_vel_afelio(), 0., 0.}); 

    GravitationalField sun_gravity(sun.get_mass(), sun.get_coordinates()); 

    SolarSystem.add_object(sun); 
    SolarSystem.add_object(earth); 

    std::cout << "Nr of planets = " << SolarSystem.get_n_objects() << std::endl; 

    for (auto i : SolarSystem.get_objects()) {
        i.print_body(); 
        i.print_position(); 
    }

    return 0; 
}
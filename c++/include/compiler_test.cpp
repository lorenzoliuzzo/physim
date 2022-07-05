// compile test

// #include "math/integral.h"
// #include "math/statistics.h"
// #include "physics/position.h"
#include "physics/system.h"
#include "physics/planet.h"
#include "physics/gravitational_field.h"


int main() {
    
    System<CelestialBody> SolarSystem; 
    Planet sun("Sun"), earth("Earth");
    Satelite moon("Moon"); 
    sun.set_position({0., 0., 0.}, {0., 0., 0.}); 
    earth.set_position({earth.get_coord_aphelion(), 0., 0.}, {earth.get_vel_aphelion(), 0., 0.}); 
    moon.set_position({earth.get_coord_aphelion() + moon.get_coord_apogee(), 0., 0.}, {earth.get_vel_aphelion() + moon.get_vel_apogee(), 0., 0.}); 
    

    GravitationalField sun_gravity(sun.get_mass(), sun.get_coordinates()); 

    SolarSystem.add_object(sun); 
    SolarSystem.add_object(earth); 
    SolarSystem.add_object(moon); 


    std::cout << "Nr of planets = " << SolarSystem.get_objects_count() << std::endl; 

    for (auto i : SolarSystem.get_objects()) {
        i.print_body(); 
    }

    return 0; 
}
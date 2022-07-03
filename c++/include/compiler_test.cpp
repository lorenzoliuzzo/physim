// compile test

// #include "math/integral.h"
// #include "math/statistics.h"
// #include "physics/position.h"
#include "physics/system.h"
#include "physics/planets.h"

int main() {
    
    System<Planet> SolarSystem; 
    Planet sun("Sun"), earth("Earth");

    SolarSystem.add_object(sun); 
    SolarSystem.add_object(earth); 

    std::cout << "Nr of planets = " << SolarSystem.get_n_objects() << std::endl; 

    for (auto i : SolarSystem.get_objects()) i.print_body(); 


    return 0; 
}
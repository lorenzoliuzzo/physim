
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Solar system.
// last updated:    05/07/2022


#include "../tools/system.h"


class SolarSystem : public system<CelestialBody> {

    public:

        void gravity(const double& h = 0.001) {
            for (unsigned int i{}; i < m_objects_count; i++) {
                GravitationalField gravity_source(m_objects[i]);    
                for (unsigned int j{}; j < m_objects_count; j++) {
                    if (i != j) {
                        appo += gravity_source.rk4(m_objects[j].get_position(), h)); 
                    }                  
                }

            }

                 

        }


};
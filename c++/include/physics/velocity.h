
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Velocity(class) for keeping track of an object that it is moving in a 3D system.
// last updated:    09/07/2022


#pragma once
#include "../math/vector_algebra.h"
#include "../physics/um.h"


class Velocity : public UM {

    protected: 

        // =============================================
        // class members
        // =============================================
    
        // velocity:     [x] [y] [z] 
        
        std::vector<double> m_velocity = zeros(3);

    public:  

        // =============================================
        // constructors
        // =============================================

        Velocity() : UM("m/s") {}

        Velocity(const std::vector<double>& vel, const char* udm_prefix = "") : m_velocity{vel}, UM("m/s", udm_prefix) {}

        
        // =============================================
        // set methods
        // =============================================

        void set_velocity(const std::vector<double>& vel) { m_velocity = vel; }

        void set_velocity_x(const double& x) { m_velocity[0] = x; }

        void set_velocity_y(const double& y) { m_velocity[1] = y;  }

        void set_velocity_z(const double& z) { m_velocity[2] = z; }
        
        
        // =============================================
        // get methods
        // =============================================

        std::vector<double> get_velocity() const { return m_velocity; }

        double get_velocity_x() const { return m_velocity[0]; }

        double get_velocity_y() const { return m_velocity[1]; }

        double get_velocity_z() const { return m_velocity[2]; }
        
        double get_magnitude() const {
            return sqrt(pow(m_velocity[0], 2) +                 
                        pow(m_velocity[1], 2) + 
                        pow(m_velocity[2], 2));
        }        

        double get_phi() const { return atan2(m_velocity[1], m_velocity[0]); }     

        double get_theta() const { return acos(m_velocity[2] / get_magnitude()); }
 
        std::vector<double> get_direction() const {
            return {cos(get_phi()), sin(get_phi()), cos(get_theta())};
        } 
       

        // =============================================
        // print methods
        // =============================================

        void print_velocity() const {
            std::cout << "velocity: ";
            UM::print_um();
            for (auto i : m_velocity) std::cout << "[" << i << "]\t";
            std::cout << std::endl; 
        }

};

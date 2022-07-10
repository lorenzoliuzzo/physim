
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Velocity(class) for keeping track of an object that it is moving in a 3D system.
// last updated:    09/07/2022


#pragma once
#include "../../math/vector_algebra.h"


class Velocity {

    protected: 

        // =============================================
        // class members
        // =============================================
    
        // velocity:     [x] [y] [z] 
        
        std::vector<double> m_velocity = zeros(3);

        const char* m_vel_um{"m"}; 

        const char* m_vel_um_prefix; 
    

    public:  

        // =============================================
        // constructors
        // =============================================

        Velocity(const char* um_prefix = "") : m_vel_um{"m/s"}, m_vel_um_prefix{um_prefix} {}

        Velocity(const std::vector<double>& vel, const char* um_prefix = "") : m_velocity{vel}, m_vel_um{"m/s"}, m_vel_um_prefix{um_prefix} {}

        
        // =============================================
        // set methods
        // =============================================

        void set_velocity(const std::vector<double>& vel) { m_velocity = vel; }

        void set_velocity_x(const double& x) { m_velocity[0] = x; }

        void set_velocity_y(const double& y) { m_velocity[1] = y;  }

        void set_velocity_z(const double& z) { m_velocity[2] = z; }
        
        void set_vel_um_prefix(const char* um_prefix) { m_vel_um_prefix = um_prefix; }


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
       
        const char* get_vel_um() const { return m_vel_um; }

        const char* get_vel_um_prefix() const { return m_vel_um_prefix; }


        // =============================================
        // print methods
        // =============================================

        void print_velocity() const {
            std::cout << "- velocity =    ";
            print(m_velocity); 
            std::cout << get_vel_um_prefix() << get_vel_um() << std::endl; 
        }

};

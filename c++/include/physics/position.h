
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Position(class) keeps track of coordinates and velocity of an object in a 3D system.
// last updated:    02/07/2022


#pragma once
#include <iostream>
#include <vector>
#include <cmath>


class Position {

    protected: 

        // =============================================
        // class members
        // =============================================
    
        // coordinates:     [i][0] 
        //    velocity:     [i][1]
        //  dimentions:   x [0][i]
        //                y [1][i] 
        //                z [2][i]

        std::vector<std::vector<double>> m_pos = zeros(3, 2); 

    
    public:  

        // =============================================
        // constructors
        // =============================================

        Position() {}

        Position(const std::vector<double>& coord, const std::vector<double>& vel) { m_pos[0] = coord; m_pos[1] = vel; }

        Position(const std::vector<std::vector<double>>& pos1) : m_pos{pos1} {}

        Position(const Position& pos1) : Position(pos1.get_position()) {}

        ~Position() {}


        // =============================================
        // set methods
        // =============================================

        void set_coordinates(const std::vector<double>& coord) { m_pos[0] = coord; }

        void set_coord_x(double x) { m_pos[0][0] = x; }

        void set_coord_y(double y) { m_pos[0][1] = y;  }

        void set_coord_z(double z) { m_pos[0][2] = z; }

        void set_velocity(const std::vector<double>& vel) { m_pos[1] = vel; }

        void set_vel_x(double vx) { m_pos[1][0] = vx; }

        void set_vel_y(double vy) { m_pos[1][1] = vy; }

        void set_vel_z(double vz) { m_pos[1][2] = vz; }

        void set_position(const std::vector<double>& coord, const std::vector<double>& vel) { m_pos[0] = coord; m_pos[1] = vel; }

        void set_position(const Position& pos1) { m_pos = pos1.get_position(); }

        void set_position(const std::vector<std::vector<double>>& pos1) { m_pos = pos1; }


        // =============================================
        // get methods
        // =============================================

        std::vector<double> get_coordinates() const { return m_pos[0]; }

        double get_coord_x() const { return m_pos[0][0]; }

        double get_coord_y() const { return m_pos[0][1]; }

        double get_coord_z() const { return m_pos[0][2]; }

        std::vector<double> get_velocity() const { return m_pos[1]; }

        double get_vel_x() const { return m_pos[1][0]; }

        double get_vel_y() const { return m_pos[1][1]; }

        double get_vel_z() const { return m_pos[1][2]; }

        std::vector<std::vector<double>> get_position() const { return m_pos; }

        double get_magnitude() const {
            return sqrt(pow(m_pos[0][0], 2) +                 
                        pow(m_pos[0][1], 2) + 
                        pow(m_pos[0][2], 2));
        }

        double get_magnitude_vel() const {
            return sqrt(pow(m_pos[1][0], 2) + 
                        pow(m_pos[1][1], 2) + 
                        pow(m_pos[1][2], 2));
        }

        double get_distance(const Position& pos1) const {
            return sqrt(pow(pos1.get_coord_x() - m_pos[0][0], 2) + 
                        pow(pos1.get_coord_y() - m_pos[0][1], 2) + 
                        pow(pos1.get_coord_z() - m_pos[0][2], 2)); 
        }

        double get_distance(const std::vector<double>& pos1) const {
            return sqrt(pow(pos1[0] - m_pos[0][0], 2) + 
                        pow(pos1[1] - m_pos[0][1], 2) + 
                        pow(pos1[2] - m_pos[0][2], 2)); 
        }

        double get_rho() { return sqrt(pow(m_pos[0][0], 2) + pow(m_pos[0][1], 2)); }

        double get_phi_coord() const { return atan2(m_pos[0][1], m_pos[0][0]); }     

        double get_phi_coord(const std::vector<double>& pos1) const { return atan2(pos1[1] - m_pos[0][1], pos1[0] - m_pos[0][0]); }
        
        double get_phi_vel() const { return atan2(m_pos[1][1], m_pos[1][0]); }

        double get_phi_vel(const std::vector<double>& vel1) const { return atan2(vel1[1] - m_pos[1][1], vel1[0] - m_pos[1][0]); }

        double get_theta_coord() const { return acos(m_pos[0][2] / get_magnitude()); }
 
        double get_theta_coord(const std::vector<double>& pos1) const { return acos((pos1[2] - m_pos[0][2]) / get_distance(pos1)); }

        double get_theta_vel() const { return acos(m_pos[1][2] / get_magnitude()); }
 
        double get_theta_vel(const std::vector<double>& vel1) const { return acos((vel1[2] - m_pos[1][2]) / get_distance(vel1)); }
 
        std::vector<double> get_coord_direction() const {
            std::vector<double> appo{cos(get_phi_coord()), sin(get_phi_coord()), cos(get_theta_coord())};
            return appo;
        } 

        std::vector<double> get_coord_direction(const std::vector<double>& coord1) const {
            std::vector<double> appo{cos(get_phi_coord(coord1)), sin(get_phi_coord(coord1)), cos(get_theta_coord(coord1))};
            return appo;
        } 

        std::vector<double> get_vel_direction() const {
            std::vector<double> appo{cos(get_phi_vel()), sin(get_phi_vel()), cos(get_theta_vel()) };
            return appo;
        } 


        std::vector<double> get_vel_direction(const std::vector<double>& vel1) const {
            std::vector<double> appo{cos(get_phi_vel(vel1)), sin(get_phi_vel(vel1)), cos(get_theta_vel(vel1))};
            return appo;
        } 

 
        // =============================================
        // print methods
        // =============================================

        void print_coordinates() const {
            std::cout << "\nCoordinates:  ";
            for (auto i : m_pos[0]) std::cout << "[" << i << "] ";
            std::cout << std::endl;
        }

        void print_velocity() const {
            std::cout << "\nVelocity:  ";
            for (auto i : m_pos[1]) std::cout << "[" << i << "] ";
            std::cout << std::endl;        
        }

        void print_position() const {
            print_coordinates(); 
            print_velocity();
        }

}; 


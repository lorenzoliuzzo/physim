
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Position(class) keeps track of coordinates and velocity of an object in a 3D system.
// last updated:    08/07/2022


#pragma once
#include "../math/vector_algebra.h"


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

        Position(const std::vector<double>& coord) { m_pos[0] = coord; }

        Position(const std::vector<double>& coord, const std::vector<double>& vel) { m_pos[0] = coord; m_pos[1] = vel; }

        Position(const std::vector<std::vector<double>>& pos1) : m_pos{pos1} {}

        Position(const Position& pos1) : Position(pos1.get_position()) {}

        ~Position() {}


        // =============================================
        // set methods
        // =============================================

        inline void set_coordinates(const std::vector<double>& coord) { m_pos[0] = coord; }

        inline void set_coord_x(double x) { m_pos[0][0] = x; }

        inline void set_coord_y(double y) { m_pos[0][1] = y;  }

        inline void set_coord_z(double z) { m_pos[0][2] = z; }

        inline void set_velocity(const std::vector<double>& vel) { m_pos[1] = vel; }

        inline void set_vel_x(double vx) { m_pos[1][0] = vx; }

        inline void set_vel_y(double vy) { m_pos[1][1] = vy; }

        inline void set_vel_z(double vz) { m_pos[1][2] = vz; }

        inline void set_position(const std::vector<double>& coord, const std::vector<double>& vel) { m_pos[0] = coord; m_pos[1] = vel; }

        inline void set_position(const Position& pos1) { m_pos = pos1.get_position(); }

        inline void set_position(const std::vector<std::vector<double>>& pos1) { m_pos = pos1; }


        // =============================================
        // get methods
        // =============================================

        inline std::vector<double> get_coordinates() const { return m_pos[0]; }

        inline double get_coord_x() const { return m_pos[0][0]; }

        inline double get_coord_y() const { return m_pos[0][1]; }

        inline double get_coord_z() const { return m_pos[0][2]; }

        inline std::vector<double> get_velocity() const { return m_pos[1]; }

        inline double get_vel_x() const { return m_pos[1][0]; }

        inline double get_vel_y() const { return m_pos[1][1]; }

        inline double get_vel_z() const { return m_pos[1][2]; }

        inline std::vector<std::vector<double>> get_position() const { return m_pos; }

        inline double get_magnitude() const {
            return sqrt(pow(m_pos[0][0], 2) +                 
                        pow(m_pos[0][1], 2) + 
                        pow(m_pos[0][2], 2));
        }

        inline double get_magnitude_vel() const {
            return sqrt(pow(m_pos[1][0], 2) + 
                        pow(m_pos[1][1], 2) + 
                        pow(m_pos[1][2], 2));
        }

        inline double get_distance(const Position& pos1) const {
            return sqrt(pow(pos1.get_coord_x() - m_pos[0][0], 2) + 
                        pow(pos1.get_coord_y() - m_pos[0][1], 2) + 
                        pow(pos1.get_coord_z() - m_pos[0][2], 2)); 
        }

        inline double get_distance(const std::vector<double>& pos1) const {
            return sqrt(pow(pos1[0] - m_pos[0][0], 2) + 
                        pow(pos1[1] - m_pos[0][1], 2) + 
                        pow(pos1[2] - m_pos[0][2], 2)); 
        }

        inline double get_rho() { return sqrt(pow(m_pos[0][0], 2) + pow(m_pos[0][1], 2)); }

        inline double get_phi_coord() const { return atan2(m_pos[0][1], m_pos[0][0]); }     

        inline double get_phi_coord(const std::vector<double>& pos1) const { return atan2(pos1[1] - m_pos[0][1], pos1[0] - m_pos[0][0]); }
        
        inline double get_phi_vel() const { return atan2(m_pos[1][1], m_pos[1][0]); }

        inline double get_phi_vel(const std::vector<double>& vel1) const { return atan2(vel1[1] - m_pos[1][1], vel1[0] - m_pos[1][0]); }

        inline double get_theta_coord() const { return acos(m_pos[0][2] / get_magnitude()); }
 
        inline double get_theta_coord(const std::vector<double>& pos1) const { return acos((pos1[2] - m_pos[0][2]) / get_distance(pos1)); }

        inline double get_theta_vel() const { return acos(m_pos[1][2] / get_magnitude()); }
 
        inline double get_theta_vel(const std::vector<double>& vel1) const { return acos((vel1[2] - m_pos[1][2]) / get_distance(vel1)); }
 
        inline std::vector<double> get_coord_direction() const {
            return {cos(get_phi_coord()), sin(get_phi_coord()), cos(get_theta_coord())};
        } 

        inline std::vector<double> get_coord_direction(const std::vector<double>& coord1) const {
            return {cos(get_phi_coord(coord1)), sin(get_phi_coord(coord1)), cos(get_theta_coord(coord1))};
        } 

        inline std::vector<double> get_vel_direction() const {
            return {cos(get_phi_vel()), sin(get_phi_vel()), cos(get_theta_vel())};
        } 


        inline std::vector<double> get_vel_direction(const std::vector<double>& vel1) const {
            return {cos(get_phi_vel(vel1)), sin(get_phi_vel(vel1)), cos(get_theta_vel(vel1))};
        } 

 
        // =============================================
        // print methods
        // =============================================

        void print_coordinates() const {
            std::cout << "Coordinates:  ";
            for (auto i : m_pos[0]) std::cout << "[" << i << "]" << std::setw(10);
            std::cout << std::endl;
        }

        void print_velocity() const {
            std::cout << "Velocity:     ";
            for (auto i : m_pos[1]) std::cout << "[" << i << "]" << std::setw(10);
            std::cout << std::endl;        
        }

        void print_position() const {
            print_coordinates(); 
            print_velocity();
        }

}; 

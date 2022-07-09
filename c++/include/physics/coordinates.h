
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Coordinates(class) for placing and keeping track of an object in a 3D system.
// last updated:    09/07/2022


#pragma once
#include "../math/vector_algebra.h"
#include "../physics/um.h"


class Coordinates : public UM {

    protected: 

        // =============================================
        // class members
        // =============================================
    
        // coordinates:     [x] [y] [z] 
        
        std::vector<double> m_coordinates = zeros(3);

    public:  

        // =============================================
        // constructors
        // =============================================

        Coordinates() : UM("m") {}

        Coordinates(const std::vector<double>& coord, const char* udm_prefix = "") : m_coordinates{coord}, UM("m", udm_prefix) {}

        
        // =============================================
        // set methods
        // =============================================

        void set_coordinates(const std::vector<double>& coord) { m_coordinates = coord; }
        
        void set_coordinate_x(const double& x) { m_coordinates[0] = x; }

        void set_coordinate_y(const double& y) { m_coordinates[1] = y;  }

        void set_coordinate_z(const double& z) { m_coordinates[2] = z; }
        
        
        // =============================================
        // get methods
        // =============================================

        std::vector<double> get_coordinates() const { return m_coordinates; }

        double get_coordinate_x() const { return m_coordinates[0]; }

        double get_coordinate_y() const { return m_coordinates[1]; }

        double get_coordinate_z() const { return m_coordinates[2]; }
        
        double get_magnitude() const {
            return sqrt(pow(m_coordinates[0], 2) +                 
                        pow(m_coordinates[1], 2) + 
                        pow(m_coordinates[2], 2));
        }        

        double get_distance(const std::vector<double>& coord) const {        
            return sqrt(pow(coord[0] - m_coordinates[0], 2) + 
                        pow(coord[1] - m_coordinates[1], 2) + 
                        pow(coord[2] - m_coordinates[2], 2)); 
        }
        
        double get_rho() { return sqrt(pow(m_coordinates[0], 2) + pow(m_coordinates[1], 2)); }

        double get_phi() const { return atan2(m_coordinates[1], m_coordinates[0]); }     

        double get_phi(const std::vector<double>& coord) const { return atan2(coord[1] - m_coordinates[1], coord[0] - m_coordinates[0]); }

        double get_theta() const { return acos(m_coordinates[2] / get_magnitude()); }
 
        double get_theta(const std::vector<double>& coord) const { return acos((coord[2] - m_coordinates[2]) / get_distance(coord)); }

        std::vector<double> get_direction() const {
            return {cos(get_phi()), sin(get_phi()), cos(get_theta())};
        } 

        std::vector<double> get_direction(const std::vector<double>& coord1) const {
            return {cos(get_phi(coord1)), sin(get_phi(coord1)), cos(get_theta(coord1))};
        } 
        

        // =============================================
        // print methods
        // =============================================

        void print_coordinates() const {
            std::cout << "coordinates: ";
            UM::print_um();
            for (auto i : m_coordinates) std::cout << "[" << i << "]\t";
            std::cout << std::endl; 
        }

};
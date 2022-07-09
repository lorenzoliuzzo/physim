
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Position(class) keeps track of coordinates and velocity of an object in a 3D system.
// last updated:    09/07/2022


#pragma once
#include "../math/vector_algebra.h"
#include "coordinates.h"
#include "velocity.h"


class Position : public Coordinates, Velocity {

    public:  

        // =============================================
        // constructors
        // =============================================

        Position() {}

        Position(const std::vector<double>& coord, const std::vector<double>& vel) : Coordinates(coord), Velocity(vel) {}

        Position(const std::vector<std::vector<double>>& pos) : Coordinates(pos[0]), Velocity(pos[1]) {}
        
        Position(const Coordinates& coord, const Velocity& vel) : Coordinates(coord), Velocity(vel) {}
        
        Position(const Position& pos) : Position(pos.get_position()) {}

        ~Position() {}


        // =============================================
        // set methods
        // =============================================

        void set_position(const std::vector<double>& coord, const std::vector<double>& vel) { 
            set_coordinates(coord); 
            set_velocity(vel); 
        }

        void set_position(const std::vector<std::vector<double>>& pos) { 
            set_coordinates(pos[0]); 
            set_velocity(pos[1]); 
        }

        void set_position(const Position& pos) { set_position(pos.get_position()); }


        // =============================================
        // get methods
        // =============================================

        std::vector<std::vector<double>> get_position() const { return {get_coordinates(), get_velocity()}; }

 
        // =============================================
        // print methods
        // =============================================

        void print_position() const {
            std::cout << "position: " << std::endl;
            print_coordinates(); 
            print_velocity();
        }

}; 

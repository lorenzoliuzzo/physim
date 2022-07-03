
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     CelestialBody(class).
// last updated:    04/07/2022


#pragma once
#include "position.h"        


class CelestialBody : public Position {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_mass, m_radius;
        const char* m_name; 


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        CelestialBody(const char* name) : Position(), m_name{name} {}

        CelestialBody(const char* name, Position pos, const double& mass, const double& radius = 0) : 
            Position(pos), m_mass{mass}, m_radius{radius}, m_name{name} {}
        
        CelestialBody(const char* name, const std::vector<double>& coord, const std::vector<double>& vel, const double& mass, const double& radius = 0) : 
            Position(coord, vel), m_mass{mass}, m_radius{radius}, m_name{name} {}
        
        CelestialBody(const char* name, const std::vector<std::vector<double>>& pos, const double& mass, const double& radius = 0) : 
            Position(pos), m_mass{mass}, m_radius{radius}, m_name{name} {}

        ~CelestialBody() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_name(const char* name) { m_name = name; }

        void set_mass(double mass) { m_mass = mass; }

        void set_radius(double r) { m_radius = r; }

        const char* get_name() const { return m_name; }

        double get_mass() const { return m_mass; }

        double get_radius() const { return m_radius; }

        
        // =============================================
        // print methods
        // =============================================

        void print_body() const {
            std::cout << "\nCelestial body: \n"; 
            std::cout << "Name = " << get_name() << std::endl; 
            std::cout << "Mass = " << get_mass() << std::endl; 
            std::cout << "Radius = " << get_radius() << std::endl;             
        }

};


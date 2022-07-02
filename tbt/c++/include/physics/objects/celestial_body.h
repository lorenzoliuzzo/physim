
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     CelestialBody(class).
// last updated:    19/06/2022

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
        
        CelestialBody() {}

        CelestialBody(const char* name, Posizione pos, const double& mass, const double& radius = 0) : 
            Posizione(pos), m_mass{mass}, m_radius{radius}, m_name{name} {}
        
        CelestialBody(const char* name, const std::vector<double>& coord, const std::vector<double>& vel, const double& mass, const double& radius = 0) : 
            Posizione(coord, vel), m_mass{mass}, m_radius{radius}, m_name{name} {}
        
        CelestialBody(const char* name, const std::vector<std::vector<double>>& pos, const double& mass, const double& radius = 0) : 
            Posizione(pos), m_mass{mass}, m_radius{radius}, m_name{name} {}

        ~CelestialBody() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_name(const char* name) { m_name = name; }

        void set_mass(double mass) { m_mass = mass; }

        void set_radius(double r) { m_radius = r; }

        void get_name() { return m_name; }

        void get_mass() { return m_mass; }

        void get_radius() { return m_radius; }

        
        // =============================================
        // print methods
        // =============================================

        void print_body() {
            std::cout << "Celestial body: \n"; 
            std::cout << "Name = " << getName() << std::endl; 
            std::cout << "Mass = " << getMass() << std::endl; 
            std::cout << "Radius = " << getRadius() << std::endl;             
            Posizione::print(); 
        }

};



// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Particle(class) is an object class that stores information about particles.
// last updated:    17/06/2022


#pragma once
#include "../tools/position.h"

#define   e_mass  9.109383701528e-31     
#define   p_mass  1.672621923695e-27
#define e_charge  -1.602176634e-19
#define p_charge  1.602176634e-19


class Particle : public Position {

    protected: 
    
        // =============================================
        // class members
        // =============================================
    
        double m_mass, m_charge;

    public:  

        // =============================================
        // constructors
        // =============================================     

        Particle(const double& mass, const double& charge, std::vector<double> coord, std::vector<double> vel) : 
            m_mass{mass}, m_charge{charge}, Position(coord, vel) {}

        Particle(const double& mass, const double& charge, Position pos) : 
            m_mass{mass}, m_charge{charge}, Position(pos) {}

        ~Particle() {}


        // =============================================
        // set methods
        // =============================================

        void set_mass(const double& mass) { m_mass = mass; }  

        void set_charge(const double& charge) { m_charge = charge; } 


        // =============================================
        // get methods
        // =============================================

        double get_mass() const { return m_mass; }

        double get_charge() const { return m_charge; }


        // =============================================
        // print methods
        // =============================================

        void print_particle() const {
            std::cout << "Mass: " << m_mass << std::endl; 
            std::cout << "Charge: " << m_charge << std::endl; 
            Position::print_position(); 
        }

};


class Electron : public Particle {

    public:  

        // =============================================
        // constructors
        // =============================================     

        Electron(std::vector<double> coord, std::vector<double> vel) : Particle(e_mass, e_charge, coord, vel) {}

        Electron(Position pos) : Particle(e_mass, e_charge, pos) {}
    
};


class Proton : public Particle {

    public:  

        // =============================================
        // constructors
        // =============================================     

        Proton(std::vector<double> coord, std::vector<double> vel) : Particle(p_mass, p_charge, coord, vel) {}

        Proton(Position pos) : Particle(p_mass, p_charge, pos) {}
    
};

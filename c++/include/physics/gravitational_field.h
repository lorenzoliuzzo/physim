
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Gravitational field.
// last updated:    05/07/2022


#pragma once
#include "ode.h"
#include "celestial_body.h"


#define G 6.6743015e-20 // udm = [km^3 kg^-1 s^-1]


class GravitationalField : public ODE {

    protected: 

        // =============================================
        // class member
        // =============================================

        double m_mass;
        Position m_pos;
        

    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        GravitationalField(const double& mass, const std::vector<double>& coord) : m_mass{mass}, m_pos(coord) {}
        
        GravitationalField(const CelestialBody& body) : m_mass{body.get_mass()}, m_pos(body.get_coordinates()) {}

        ~GravitationalField() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_mass(const double& mass) { m_mass = mass; }

        double get_mass() const { return m_mass; }

        void set_coord(const std::vector<double>& coord) { m_pos.set_coordinates(coord); }

        std::vector<double> get_coord() const { return m_pos.get_coordinates(); }
        

        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) override {
            double phi{m_pos.get_phi_coord(pos[0])}; 
            std::vector<double> direction{cos(phi), sin(phi), 0}; 
            m_df[0] = pos[1]; 
            m_df[1] = direction * (- G * m_mass / pow(m_pos.get_distance(pos[0]), 2));
            return m_df; 
        }

};


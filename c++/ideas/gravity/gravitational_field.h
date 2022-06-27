
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     GravitationalField(class) is an ode solver based on Newton's general equations of gravitation.
// last updated:    19/06/2022

#pragma once
#include "ode.h"
#include "system.h"

#define G 6.6743015e-11


class GravitationalField : public ODE {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_mass;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        GravitationalField(const SystemBase& system) : 
            ODE(coord, vel, t), m_mass{mass} {}

        GravitationalField(const std::vector<std::vector<double>>& pos, double mass, double t = 0) : 
            ODE(pos, t), m_mass{mass} {}

        GravitationalField(const Position& pos, double t = 0) : 
            ODE(pos, t), m_mass{mass} {}
        
        ~GravitationalField() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_mass1(double mass1) { m_mass1 = mass1; }

        void set_mass2(double mass2) { m_mass2 = mass2; }

        double get_mass1() const { return m_mass1; }

        double get_mass2() const { return m_mass2; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& init, double h = 0.001) override {
            reset_df(); 
            m_df[0] = init[1]; 
            m_df[1] -= G * ; 
            return m_df;
        }


};}


#pragma once
#include "../../eq_differenziali.h"
#include "../../objects/celestial_body.h"

#define G 6.6742e-20

class CampoGravitazionale : public FunzioneVettorialeBase, Posizione {
 private: 
    double m_massa; 

 public: 
    CampoGravitazionale() : m_massa{}, Posizione() {}

    CampoGravitazionale(double massa, Coordinate coord) : 
        m_massa{massa}, Posizione(coord) {}

    CampoGravitazionale(double massa, std::vector<double> x) : 
        m_massa{massa}, Posizione(x) {}

    CampoGravitazionale(CelestialBody p) :
        m_massa{p.getMass()}, 
        Posizione(p.getPos()) {}

    // CampoGravitazionale(double massa, double x, double y, double z) : 
    //     m_massa{massa}, Posizione(x, y, z) {}

    ~CampoGravitazionale() {}

    void setMassa(double m) { m_massa = m; } 
    double getMassa() const { return m_massa; }

    Matrix eval(const Matrix & x, double t) const override {
        Matrix result(x.getRowSize(), x.getColSize()); 
        double phi{getPhi(x.getCol(0))}; 
        double pippo = - G * m_massa / pow(getDistanza(x.getCol(0)), 2); 
        std::vector<double> appo{cos(phi), sin(phi), 0}; 
        appo = appo * pippo; 
        result.setCol(x.getCol(1), 0);
        result.setCol(appo, 1); 
        return result; 
    }
};

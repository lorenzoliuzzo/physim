
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Different types of oscillators.
// last updated:    19/06/2022

#pragma once
#include "../ode.h"


class HarmonicOscillator : public ODE {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_omega;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        HarmonicOscillator(const std::vector<double>& coord, const std::vector<double>& vel, double omega, double t = 0) : 
            ODE(coord, vel, t), m_omega{omega} {}

        HarmonicOscillator(const std::vector<std::vector<double>>& pos, double omega, double t = 0) : 
            ODE(pos, t), m_omega{omega} {}

        HarmonicOscillator(const Position& pos, double omega, double t = 0) : 
            ODE(pos, t), m_omega{omega} {}
        
        ~HarmonicOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega(double omega) { m_omega = omega; }

        double get_omega() const { return m_omega; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& init, double h = 0.001) override {
            reset_df(); 
            m_df[0] = init[1]; 
            m_df[1] = - pow(m_omega, 2) * init[0]; 
            return m_df;
        }


};

class ForcedOscillator : public ODE {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_omega0, m_omega1;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        ForcedOscillator(const std::vector<double>& coord, const std::vector<double>& vel, double omega0, double omega1, double t = 0) : 
            ODE(coord, vel, t), m_omega0{omega0}, m_omega1{omega1} {}

        ForcedOscillator(const std::vector<std::vector<double>>& pos, double omega0, double omega1, double t = 0) : 
            ODE(pos, t), m_omega0{omega0}, m_omega1{omega1} {}

        ForcedOscillator(const Position& pos, double omega0, double omega1, double t = 0) : 
            ODE(pos, t), m_omega0{omega0}, m_omega1{omega1} {}
        
        ~ForcedOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega0(double omega0) { m_omega0 = omega0; }

        void set_omega1(double omega1) { m_omega1 = omega1; }

        double get_omega0() const { return m_omega0; }

        double get_omega1() const { return m_omega1; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& init, double h = 0.001) override {
            reset_df(); 
            m_df[0] = init[1]; 
            m_df[1] += (- pow(m_omega0, 2) * init[0] + sin(m_omega1 * m_time)); 
            return m_df; 
        }


};


class DampedOscillator : public ODE {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_omega, m_alpha;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        DampedOscillator(const std::vector<double>& coord, const std::vector<double>& vel, double omega, double alpha, double t = 0) : 
            ODE(coord, vel, t), m_omega{omega}, m_alpha{alpha} {}

        DampedOscillator(const std::vector<std::vector<double>>& pos, double omega, double alpha, double t = 0) : 
            ODE(pos, t), m_omega{omega}, m_alpha{alpha} {}

        DampedOscillator(const Position& pos, double omega, double alpha, double t = 0) : 
            ODE(pos, t), m_omega{omega}, m_alpha{alpha} {}
        

        ~DampedOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega(double omega) { m_omega = omega; }

        void set_alpha(double alpha) { m_alpha = alpha; }

        double get_omega() const { return m_omega; }

        double get_alpha() const { return m_alpha; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& init, double h = 0.001) override {
            reset_df(); 
            m_df[0] = init[1]; 
            m_df[1] -= ((pow(m_omega, 2) * init[0]) + (m_alpha * init[1])); 
            return m_df; 
        }


};


class ForcedDampedOscillator : public ODE {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_omega0, m_omega1, m_alpha;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        ForcedDampedOscillator(const std::vector<double>& coord, const std::vector<double>& vel, double omega0, double omega1, double alpha, double t = 0) : 
            ODE(coord, vel, t), m_omega0{omega0}, m_omega1{omega1}, m_alpha{alpha} {}

        ForcedDampedOscillator(const std::vector<std::vector<double>>& pos, double omega0, double omega1, double alpha, double t = 0) : 
            ODE(pos, t), m_omega0{omega0}, m_omega1{omega1}, m_alpha{alpha} {}

        ForcedDampedOscillator(const Position& pos, double omega0, double omega1, double alpha, double t = 0) : 
            ODE(pos, t), m_omega0{omega0}, m_omega1{omega1}, m_alpha{alpha} {}
        

        ~ForcedDampedOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega0(double omega0) { m_omega0 = omega0; }

        void set_omega1(double omega1) { m_omega1 = omega1; }

        void set_alpha(double alpha) { m_alpha = alpha; }

        double get_omega0() const { return m_omega0; }

        double get_omega1() const { return m_omega1; }

        double get_alpha() const { return m_alpha; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& init, double h = 0.001) override {
            reset_df(); 
            m_df[0] = init[1]; 
            m_df[1] += (sin(m_omega1 * m_time) - pow(m_omega0, 2) * init[0] - m_alpha * init[1]); 
            return m_df; 
        }


};
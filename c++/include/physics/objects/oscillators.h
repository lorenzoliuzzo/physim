
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Different types of oscillators.
// last updated:    02/07/2022


#pragma once
#include "../../math/ode.h"


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
        
        HarmonicOscillator(const double& omega, const double& time = 0) : m_omega{omega} { Time::m_time = time; }
        
        ~HarmonicOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega(const double& omega) { m_omega = omega; }

        double get_omega() const { return m_omega; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) override {
            m_df[0] = pos[1]; 
            m_df[1] = - pow(m_omega, 2) * pos[0]; 
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
        
        ForcedOscillator(const double& omega0, const double& omega1, const double& time = 0) : 
            m_omega0{omega0}, m_omega1{omega1} { Time::m_time = time; }
        
        ~ForcedOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega0(const double& omega0) { m_omega0 = omega0; }

        void set_omega1(const double& omega1) { m_omega1 = omega1; }

        double get_omega0() const { return m_omega0; }

        double get_omega1() const { return m_omega1; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) override {
            m_df[0] = pos[1]; 
            m_df[1] = - pow(m_omega0, 2) * pos[0] + sin(m_omega1 * get_time()); 
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
        
        DampedOscillator(const double& omega, const double& alpha, const double& time = 0) : 
            m_omega{omega}, m_alpha{alpha} { Time::m_time = time; }

        ~DampedOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega(const double& omega) { m_omega = omega; }

        void set_alpha(const double& alpha) { m_alpha = alpha; }

        double get_omega() const { return m_omega; }

        double get_alpha() const { return m_alpha; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) override {
            m_df[0] = pos[1]; 
            m_df[1] = - pow(m_omega, 2) * pos[0] - m_alpha * pos[1]; 
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
        
        ForcedDampedOscillator(const double& omega0, const double& omega1, const double& alpha, const double& time = 0) : 
            m_omega0{omega0}, m_omega1{omega1}, m_alpha{alpha} { Time::m_time = time; }

        ~ForcedDampedOscillator() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_omega0(const double& omega0) { m_omega0 = omega0; }

        void set_omega1(const double& omega1) { m_omega1 = omega1; }

        void set_alpha(const double& alpha) { m_alpha = alpha; }

        double get_omega0() const { return m_omega0; }

        double get_omega1() const { return m_omega1; }

        double get_alpha() const { return m_alpha; }


        // =============================================
        // eval methods
        // =============================================

        std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) override {
            m_df[0] = pos[1]; 
            m_df[1] = - pow(m_omega0, 2) * pos[0] - m_alpha * pos[1] + sin(m_omega1 * get_time()); 
            return m_df; 
        }


};

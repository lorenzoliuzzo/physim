
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Ordinaries Differential Equations.
// last updated:    02/07/2022


#pragma once
#include "../math/vector_algebra.h"
#include "position.h"


class ODE {

    public: 

        // =============================================
        // class members
        // =============================================

        std::vector<std::vector<double>> m_df = zeros(3, 2); 
        
        double m_time{}, m_h; 


        // =============================================
        // virtual destructor
        // =============================================

        virtual ~ODE() {}


        // =============================================
        // time methods
        // =============================================
        
        double get_time() const { return m_time; }

        void increase_time(const double& h) { m_time += h; } 

        void reset_time() { m_time = 0; }


        // =============================================
        // virtual eval methods
        // =============================================

        virtual std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) = 0; 


        // =============================================
        // integration methods
        // =============================================

        std::vector<std::vector<double>> euler(const std::vector<std::vector<double>>& pos, const double& h = 0.001) {
            std::vector<std::vector<double>> appo = pos + h * eval(pos, m_time); 
            return appo; 
        }

        std::vector<std::vector<double>> euler_modified(const std::vector<std::vector<double>>& pos, const double& h = 0.001) {
            std::vector<std::vector<double>> appo = pos + h * eval(pos, m_time); 
            appo = pos + h * (eval(pos, m_time) + eval(appo, m_time + h)) / 2.; 
            return appo; 
        }

        std::vector<std::vector<double>> rk4(const std::vector<std::vector<double>>& pos, const double& h = 0.001) {
            std::vector<std::vector<double>> k1{}, k2{}, k3{}, k4{}; 
            k1 = eval(pos, m_time); 
            k2 = eval(pos + k1 * h / 2., m_time + h / 2.);
            k3 = eval(pos + k2 * h / 2., m_time + h / 2.);
            k4 = eval(pos + k3 * h / 2., m_time + h / 2.);      
            return (pos + (k1 + k2 * 2. + k3 * 2. + k4) * (h / 6.)); 
        } 

};


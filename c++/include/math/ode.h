
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Ordinaries Differential Equations.
// last updated:    10/07/2022


#pragma once
#include "vector_algebra.h" 
#include "../physics/tools/time.h"


class ODE : public Time {

    public: 

        // =============================================
        // class members
        // =============================================

        std::vector<std::vector<double>> m_df = zeros(3, 2); 
        
        double m_h; 


        // =============================================
        // virtual destructor
        // =============================================

        virtual ~ODE() {}


        // =============================================
        // virtual eval methods
        // =============================================

        virtual std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& pos, const double& h = 0.001) = 0; 


        // =============================================
        // integration methods
        // =============================================

        std::vector<std::vector<double>> euler(const std::vector<std::vector<double>>& pos, const double& h = 0.001) {
            std::vector<std::vector<double>> appo = pos + h * eval(pos, get_time()); 
            return appo; 
        }

        std::vector<std::vector<double>> euler_modified(const std::vector<std::vector<double>>& pos, const double& h = 0.001) {
            std::vector<std::vector<double>> appo = pos + h * eval(pos, get_time()); 
            appo = pos + h * (eval(pos, get_time()) + eval(appo, get_time() + h)) / 2.; 
            return appo; 
        }

        std::vector<std::vector<double>> rk4(const std::vector<std::vector<double>>& pos, const double& h = 0.001) {
            std::vector<std::vector<double>> k1{}, k2{}, k3{}, k4{}; 
            k1 = eval(pos, get_time()); 
            k2 = eval(pos + k1 * h / 2., get_time() + h / 2.);
            k3 = eval(pos + k2 * h / 2., get_time() + h / 2.);
            k4 = eval(pos + k3 * h / 2., get_time() + h / 2.);      
            return (pos + (k1 + k2 * 2. + k3 * 2. + k4) * (h / 6.)); 
        } 

};



// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Ordinaries Differential Equations.
// last updated:    19/06/2022

#pragma once
#include "../math/vector_algebra.h"
#include "position.h"
#include <chrono>


class ODE : public Position {

    protected: 

        // =============================================
        // class members
        // =============================================

        std::vector<std::vector<double>> m_df; 
        double m_time, m_h; 


    public: 

        // =============================================
        // constructor 
        // =============================================

        ODE(const std::vector<double>& coord, const std::vector<double>& vel, double t = 0) : Position(coord, vel), m_time{t} {
            m_df.resize(2);
            for (int i{}; i < 3; i++) {
                m_df[i].resize(3);
            }
        }

        ODE(const std::vector<std::vector<double>>& pos, double t = 0) : Position(pos), m_time{t} {
            m_df.resize(2);
            for (int i{}; i < 3; i++) {
                m_df[i].resize(3);
            }
        }

        ODE(const Position& pos, double t = 0) : Position(pos), m_time{t} {
            m_df.resize(2);
            for (int i{}; i < 3; i++) {
                m_df[i].resize(3);
            }
        }
        ~ODE() {}


        // =============================================
        // useful methods
        // =============================================
        
        double get_time() const { return m_time; }

        void reset_time() { m_time = 0; }

        void reset_df() { 
            for (size_t i = 0; i < m_df.size(); i++) {
                for (size_t j = 0; j < m_df[i].size(); j++) {
                    m_df[i][j] = 0;
                }
            }
        }


        // =============================================
        // eval methods
        // =============================================

        virtual std::vector<std::vector<double>> eval(const std::vector<std::vector<double>>& init, double h = 0.001) = 0; 


        // =============================================
        // integration methods
        // =============================================

        void euler(double h = 0.001) {
            std::vector<std::vector<double>> appo = get_position();
            set_position(appo + h * eval(appo, m_time)); 
            m_time += h; 
        }

        void rk4(double h = 0.001) {
            std::vector<std::vector<double>> k1{}, k2{}, k3{}, k4{}; 
            k1 = eval(get_position(), m_time); 
            k2 = eval(get_position() + k1 * h / 2., m_time + h / 2.);
            k3 = eval(get_position() + k2 * h / 2., m_time + h / 2.);
            k4 = eval(get_position() + k3 * h / 2., m_time + h / 2.);             
            set_position(get_position() + (k1 + k2 * 2. + k3 * 2. + k4) * (h / 6.)); 
            m_time += h; 
        } 


};


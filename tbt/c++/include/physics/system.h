
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     GravitationalField(class) is an ode solver based on Newton's general equations of gravitation.
// last updated:    19/06/2022

#pragma once
#include "ode.h"


template <typename T> 
class SystemBase {

    private: 

        // =============================================
        // class member
        // =============================================

        std::vector<T> m_objects; 
        unsigned int m_objects_count; 
        double m_kinetic_energy, m_potential_energy, m_total_energy; 


    public:

        // =============================================
        // Constructor and destructor
        // =============================================

        SystemBase() :
            m_kinetic_energy{},
            m_potential_energy{}, 
            m_total_energy{} {} 

        ~SystemBase() {}


        // =============================================
        // add and delete methods
        // =============================================
     
        void add_object(T x) { m_objects.push_back(x); } 

        void del_object(unsigned int pos) { m_objects.erase(pos); }


        // =============================================
        // set and get methods
        // =============================================

        void set_kinetic_energy(double K) { m_kinetic_energy = K; }

        void set_potential_energy(double U) { m_potential_energy = U; }

        void set_total_energy(double E) { m_total_energy = E; }

        void reset_objects() { m_objects.clear(); }

        void resetEnergy() {
            m_kineticEnergy = 0; 
            m_potentialEnergy = 0; 
            m_totalEnergy = 0;            
        }

        int get_objects_count() const { return m_objects.size(); }

        double get_total_energy() const { return m_total_energy; }

        double get_potential_energy() const { return m_potential_energy; }

        double get_kinetic_energy() const { return m_kinetic_energy; }

        std::vector<T> &get_objects() { return m_objects; }

        T &get_object(unsigned pos) { return m_objects[pos]; }


};

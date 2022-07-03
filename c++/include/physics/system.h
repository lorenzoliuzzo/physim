
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Basic physics system of generic objects.
// last updated:    04/07/2022


#include "../math/vector_algebra.h"
#include "position.h"


template <typename T> 
class System : Position {

    private:
        
        // =============================================
        // class members
        // =============================================     

        std::vector<T> m_objects; 
        unsigned int m_objects_count; 

        double m_kinetic_energy;
        double m_potential_energy;
        double m_total_energy; 


    public:

        // =============================================
        // constructor and destructor
        // =============================================     

        System() :
            m_kinetic_energy{},
            m_potential_energy{}, 
            m_total_energy{} {} 

        ~System() {}

        void add_object(T x) { m_objects.push_back(x); } 

        void set_kinetic_energy(double K) { m_kinetic_energy = K; }

        void set_potential_energy(double U) { m_potential_energy = U; }

        void set_total_energy(double E) { m_total_energy = E; }

        void reset_Objects() {
            m_objects.clear(); 
        }

        void reset_Energy() {
            m_kinetic_energy = 0; 
            m_potential_energy = 0; 
            m_total_energy = 0;            
        }

        int get_n_objects() const {
            return m_objects.size();
        }

        double get_total_energy() const {
            return m_total_energy;
        }

        double get_potential_energy() const {
            return m_potential_energy;
        }

        double get_kinetic_energy() const {
            return m_kinetic_energy;
        }

        std::vector<T> &get_objects() {
            return m_objects;
        }

        T &get_object(unsigned pos) {
            return m_objects[pos]; 
        }

};

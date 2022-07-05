
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Basic physics system of generic objects.
// last updated:    05/07/2022


#include "../math/vector_algebra.h"
#include "gravitational_field.h"

template <typename T> 
class System {

    private:
        
        // =============================================
        // class members
        // =============================================     

        std::vector<T> m_objects; 
        unsigned int m_objects_count{}, m_field_count{}; 
        double m_time; 
        bool m_gravitational_field, m_electric_field, m_magnetic_field;
        std::vector<GravitationalField> m_gravity; 


    public:

        // =============================================
        // constructor and destructor
        // =============================================     

        System() {}

        ~System() {}


        // =============================================
        // set & get methods
        // =============================================     
        
        void add_object(T x) { m_objects.push_back(x); } 

        void reset_objects() {
            m_objects.clear(); 
        }

        int get_objects_count() const {
            return m_objects.size();
        }

        std::vector<T> &get_objects() {
            return m_objects;
        }

        T &get_object(unsigned pos) {
            return m_objects[pos]; 
        }
        

        // =============================================
        // time methods
        // =============================================
        
        double get_time() const { return m_time; }

        void increase_time(const double& h) { m_time += h; } 

        void reset_time() { m_time = 0; }
        

        // =============================================
        // fields methods
        // =============================================
        
        void activate_gravitational_field() {
            m_gravitational_field = true; 
            m_field_count++; 
            for (unsigned int i{}; i < get_objects_count(); i++) {
                m_gravity.push_back(GravitationalField(get_object(i))); 
            }
        }        
        
        // void activate_electric_field() {
        //     m_electric_field = true; 
        //     m_field_count++; 
        //     for (unsigned int i{}; i < get_objects_count(); i++) {
        //         m_gravity.push_back(new GravitationalField(get_object(i))); 
        //     }
        // }        
        
        // void activate_magnetic_field() {
        //     m_magnetic_field = true; 
        //     m_field_count++; 
        //     for (unsigned int i{}; i < get_objects_count(); i++) {
        //         m_gravity.push_back(new GravitationalField(get_object(i))); 
        //     }
        // }


        // =============================================
        // evolve methods
        // =============================================
        
        void evolve(const double& h = 1.e-3) {
            std::vector<std::vector<double>> appo = zeros(3, 2); 
            unsigned int i{}; 
            switch(m_field_count) {
                case 0: 
                    std::cout << "There are no fields activated. Retry!" << std::endl;
                    exit(-11); 
                
                case 1: 
                    if (m_gravitational_field == true) {
                        for (auto k : m_objects) {
                            for (unsigned int j{}; j < get_objects_count(); j++) {
                                if (j != i) appo += m_gravity[j].rk4(k.get_position(), h); 
                            }
                            k.set_position(appo); 
                            increase_time(h); 
                            i++; 
                            appo.clear(); 
                        }
                    }
                // case 2:
            }
        }

        
        // void evolve_and_plot(const char* filename, const char* axis_x = "", const char* axis_y = "", const double& h = 1.e-3) {
        //     plot.redirect_to_png(filename);
        //     plot.set_xlabel(axis_x);
        //     plot.set_ylabel(axis_y);  
        //     std::vector<double> coord_x{}, coord_y{};
        //     std::vector<std::vector<double>> appo = zeros(3, 2); 

        //     unsigned int i{}; 
        //     switch(m_field_count) {
        //         case 0: 
        //             std::cout << "There are no fields activated. Retry!" << std::endl;
        //             exit(-11); 
                
        //         case 1: 
        //             if (m_gravitational_field == true) {
        //                 for (auto k : m_objects) {
        //                     for (unsigned int j{}; j < get_objects_count(); j++) {
        //                         if (j != i) appo += m_gravity[j].rk4(k.get_position(), h); 
        //                     }
        //                     k.set_position(appo); 
        //                     increase_time(h); 
        //                     i++; 
        //                     appo.clear(); 
        //                 }
        //             }
        //         // case 2:
        //     }
        // }
        
};


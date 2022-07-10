
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Mass(class) defining one of the most common object member and 
//                  the source of the gravitational field.
// last updated:    10/07/2022


#include "coordinates.h"


#define G 6.6743015e-11 // um = [m^3 kg^-1 s^-1]

 
class Mass : public Coordinates {

    protected: 

        // =============================================
        // class member
        // =============================================

        double m_mass;

        const char* m_mass_um{"g"}; 

        const char* m_mass_um_prefix; 
        
        bool m_gravitational_field = false;

        std::vector<double> m_gravitational_attraction = zeros(3); 


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        Mass(const double& mass, const char* um_prefix = "") : m_mass{mass}, m_mass_um_prefix{um_prefix}, Coordinates(zeros(3)) {}

        Mass(const double& mass, const std::vector<double>& coord, const char* mass_um_prefix = "", const char* coord_udm_prefix = "") : m_mass{mass}, m_mass_um_prefix{mass_um_prefix}, Coordinates(coord, coord_udm_prefix) {}

        Mass(const double& mass, const Coordinates& coord, const char* mass_um_prefix = "") : m_mass{mass}, m_mass_um_prefix{mass_um_prefix}, Coordinates(coord) {}

        ~Mass() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_mass(const double& mass) { m_mass = mass; }

        void set_mass_um_prefix(const char* um_prefix) { m_mass_um_prefix = um_prefix; }
        
        double get_mass() const { return m_mass; }

        const char* get_mass_um() const { return m_mass_um; }

        const char* get_mass_um_prefix() const { return m_mass_um_prefix; }

      
        // =============================================
        // gravitational methods
        // =============================================

        void activate_gravitational_field() { m_gravitational_field = true; }

        void deactivate_gravitational_field() { m_gravitational_field = false; }

        void reset_gravitational_attraction() { m_gravitational_attraction.clear(); }

        void add_gravitational_attraction(const std::vector<double>& attraction) { m_gravitational_attraction += attraction; }

        std::vector<double> get_gravitational_attraction() const { return m_gravitational_attraction; }

        std::vector<double> gravitational_attraction(const std::vector<double>& coord1) {
            if (m_gravitational_field == false) {
                std::cout << "Before evaluating the gravitational attraction given by this mass in these coordinates, you must activate the gravitational field." << std::endl; 
                exit(-11);
            }
            if (coord1 == get_coordinates()) return zeros(3);
            std::vector<double> appo = get_direction(coord1) * (- G * m_mass / pow(get_distance(coord1), 2));
            return appo; 
        }


        // =============================================
        // print methods
        // =============================================
        
        void print_mass() const { 
            std::cout << "- mass = " << get_mass() << " " << get_mass_um_prefix() << get_mass_um() << std::endl;
        }

        void print_gravitational_attraction() const { 
            std::cout << "- gravitational attraction: " << std::endl; 
            for (auto i : get_gravitational_attraction()) std::cout << "[" << i << "]\t";
            std::cout << std::endl; 
        }
    
};


// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Mass(class) defining one of the most common object member and 
//                  the source of the gravitational field.
// last updated:    08/07/2022


#include "position.h"


#define G 6.6743015e-11 // udm = [m^3 kg^-1 s^-1]


class Mass : public Position {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_mass;

        bool m_gravitational_field;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        Mass(const double& mass, const std::vector<double>& coord = zeros(3)) : m_mass{mass}, Position(coord) {}

        Mass(const double& mass, const std::vector<std::vector<double>>& pos = zeros(3, 2)) : m_mass{mass}, Position(pos) {}

        ~Mass() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_mass(const double& mass) { m_mass = mass; }

        inline double get_mass() const { return m_mass; }

        inline void print_mass() const { std::cout << "Mass = " << get_mass() << std::endl; }


        // =============================================
        // gravitational methods
        // =============================================

        inline void activate_gravitational_field() { m_gravitational_field = true; }

        inline void deactivate_gravitational_field() { m_gravitational_field = false; }

        std::vector<double> gravitational_attraction(const std::vector<double>& coord1, const char* udm = "m") {
            if (m_gravitational_field == false) {
                std::cout << "Before evaluating the gravitational attraction given by this mass in these coordinates, you must activate the gravitational field." << std::endl; 
                exit(-11);
            }
            if (coord1 == get_coordinates()) return zeros(3);
            double phi{get_phi_coord(coord1)}; 
            std::vector<double> direction{cos(phi), sin(phi), 0}; 
            std::vector<double> appo = direction * (- G * m_mass / pow(get_distance(coord1), 2));
            if (udm == "km") return appo * 1.e-9;
            if (udm == "dm") return appo * 1.e3; 
            if (udm == "cm") return appo * 1.e6; 
            if (udm == "mm") return appo * 1.e9; 
            if (udm == "microm") return appo * 1.e18;
            if (udm == "nm") return appo * 1.e27;
            else return appo;
        }

};

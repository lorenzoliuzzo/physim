
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Charge(class) defining one of the most common object member and 
//                  the source of the electric field.
// last updated:    08/07/2022


#include "position.h"


#define vacuun_permittivity 8.854187812813e-12 // udm [N^-1 m^-2 C^2]
#define K 1. / (4 * M_PI * vacuun_permittivity) // udm [N m^2 C^-2]


class Charge : public Position {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_charge;

        double m_permittivity;

        bool m_electric_field;


    public: 

        // =============================================
        // constructors and destructor
        // =============================================

        Charge(const double& charge, const std::vector<double>& coord = zeros(3), const double& permittivity = 1.) : m_charge{charge}, Position(coord), m_permittivity{permittivity} {}

        Charge(const double& charge, const std::vector<std::vector<double>>& pos = zeros(3, 2), const double& permittivity = 1.) : m_charge{charge}, Position(pos), m_permittivity{permittivity} {}

        ~Charge() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_charge(const double& charge) { m_charge = charge; }

        double get_charge() const { return m_charge; }

        void set_permittivity(const double& permittivity) { m_permittivity = permittivity; }

        double get_permittivity() const { return m_permittivity; }
        

        // =============================================
        // electric methods
        // =============================================

        void activate_electric_field() { m_electric_field = true; }

        void deactivate_electric_field() { m_electric_field = false; }

        std::vector<double> electric_attraction(const std::vector<double>& coord1, const char* udm = "C") {
            if (m_electric_field == false) {
                std::cout << "Before evaluating the electric attraction given by this charge in these coordinates, you must activate the electric field." << std::endl; 
                exit(-11);
            }
            if (coord1 == get_coordinates()) return zeros(3);
            double phi{get_phi_coord(coord1)}; 
            std::vector<double> direction{cos(phi), sin(phi), 0}; 
            std::vector<double> appo = direction * (K * m_charge / (m_permittivity * pow(get_distance(coord1), 2)));
            if (udm == "mC") return appo * 1.e6; 
            if (udm == "microC") return appo * 1.e12;
            if (udm == "nC") return appo * 1.e36;
            else return appo;
        }

};

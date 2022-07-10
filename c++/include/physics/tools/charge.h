
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Charge(class) defining one of the most common object member and 
//                  the source of the electric field.
// last updated:    10/07/2022


#include "coordinates.h"


#define vacuun_permittivity 8.854187812813e-12 // um [N^-1 m^-2 C^2]
#define K 1. / (4 * M_PI * vacuun_permittivity) // um [N m^2 C^-2]

 
class Charge : public Coordinates {

    private: 

        // =============================================
        // class member
        // =============================================

        double m_charge;

        double m_permittivity; 

        const char* m_charge_um{"C"}; 

        const char* m_charge_um_prefix; 

        bool m_elettric_field = false;

        std::vector<double> m_elettric_attraction = zeros(3); 


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        Charge(const double& charge, const double& permittivity = 1, const char* um_prefix = "") : m_charge{charge}, m_permittivity{permittivity}, m_charge_um_prefix{um_prefix}, Coordinates(zeros(3)) {}

        Charge(const double& charge, const std::vector<double>& coord, const double& permittivity = 1, const char* charge_um_prefix = "", const char* coord_udm_prefix = "") : m_charge{charge}, m_permittivity{permittivity}, m_charge_um_prefix{charge_um_prefix}, Coordinates(coord, coord_udm_prefix) {}

        Charge(const double& charge, const Coordinates& coord, const double& permittivity = 1, const char* charge_um_prefix = "") : m_charge{charge}, m_permittivity{permittivity}, m_charge_um_prefix{charge_um_prefix}, Coordinates(coord) {}

        ~Charge() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_charge(const double& charge) { m_charge = charge; }

        void set_permittivity(const double& permittivity) { m_permittivity = permittivity; }

        void set_charge_um_prefix(const char* um_prefix) { m_charge_um_prefix = um_prefix; }

        double get_charge() const { return m_charge; }

        double get_permittivity() const { return m_permittivity; }

        const char* get_charge_um() const { return m_charge_um; }

        const char* get_charge_um_prefix() const { return m_charge_um_prefix; }

        
        // =============================================
        // elettric methods
        // =============================================

        void activate_elettric_field() { m_elettric_field = true; }

        void deactivate_elettric_field() { m_elettric_field = false; }
        
        void reset_elettric_attraction() { m_elettric_attraction.clear(); }

        void set_elettric_attraction(const std::vector<double>& attraction) { m_elettric_attraction = attraction; }

        std::vector<double> get_elettric_attraction() const { return m_elettric_attraction; }

        std::vector<double> elettric_attraction(const std::vector<double>& coord1) {
            if (m_elettric_field == false) {
                std::cout << "Before evaluating the elettric attraction given by this charge in these coordinates, you must activate the elettric field." << std::endl; 
                exit(-11);
            }
            if (coord1 == get_coordinates()) return zeros(3);
            else return get_direction(coord1) * K * m_charge / (m_permittivity *  pow(get_distance(coord1), 2));
        }


        // =============================================
        // print methods
        // =============================================
        
        void print_charge() const { 
            std::cout << "charge = " << get_charge() << " " << get_charge_um_prefix() << get_charge_um() << std::endl;
        }

        void print_elettric_attraction() const { 
            std::cout << "elettric attraction: " << std::endl; 
            for (auto i : get_elettric_attraction()) std::cout << "[" << i << "]\t";
            std::cout << std::endl; 
        }

};


// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     CelestialBody(class) mother class of Planet(class), Satelites(class) and DwarfPlanet(class).
// last updated:    05/07/2022


#pragma once
#include "../tools/mass.h"   
#include "../tools/position.h"   
#include "../tools/shape.h"     


class CelestialBody : public Mass, public Position, public Sphere {

    protected: 

        // =============================================
        // class member
        // =============================================

        const char* m_name; 
        
        const char* m_type; 


    public: 

        // =============================================
        // constructors and destructor
        // =============================================
        
        CelestialBody(const char* name, const char* type) : m_name{name}, m_type{type}, Mass(0.), Position(), Sphere(0., zeros(3)) {}

        CelestialBody(const char* name, 
                    const char* type,
                    const double& mass,
                    const double& radius,
                    const std::vector<double>& coord,
                    const std::vector<double>& vel,
                    const char* mass_um_prefix = "",
                    const char* coord_um_prefix = "",
                    const char* vel_um_prefix = "", 
                    const char* shape_um_prefix = "") : 
            m_name{name}, m_type{type}, Mass(mass, mass_um_prefix), Position(coord, vel, coord_um_prefix, vel_um_prefix), Sphere(radius, coord, shape_um_prefix) {}
        
        CelestialBody(const char* name, 
                    const char* type,
                    const double& mass,
                    const double& radius,
                    const Coordinates& coord,
                    const Velocity& vel,
                    const char* mass_um_prefix = "",
                    const char* shape_um_prefix = "") : 
            m_name{name}, m_type{type}, Mass(mass, mass_um_prefix), Position(coord, vel), Sphere(radius, coord.get_coordinates(), shape_um_prefix) {}

        CelestialBody(const char* name, 
                    const char* type,
                    const double& mass,
                    const double& radius,
                    const std::vector<std::vector<double>>& pos,
                    const char* mass_um_prefix = "",
                    const char* coord_um_prefix = "",
                    const char* vel_um_prefix = "", 
                    const char* shape_um_prefix = "") : 
            m_name{name}, m_type{type}, Mass(mass, mass_um_prefix), Position(pos, coord_um_prefix, vel_um_prefix), Sphere(radius, pos[0], shape_um_prefix) {}

        CelestialBody(const char* name, 
                    const char* type,
                    const double& mass,
                    const double& radius,
                    const Position& pos,
                    const char* mass_um_prefix = "",
                    const char* shape_um_prefix = "") : 
            m_name{name}, m_type{type}, Mass(mass, mass_um_prefix), Position(pos), Sphere(radius, pos.get_coordinates(), shape_um_prefix) {}

        ~CelestialBody() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_name(const char* name) { m_name = name; }

        const char* get_name() const { return m_name; }

        const char* get_type() const { return m_type; }

        
        // =============================================
        // print methods
        // =============================================

        void print_name() const { std::cout << "- name = " << get_name() << std::endl; }
        
        void print_type() const { std::cout << "- type = " << get_type() << std::endl; }

        void print_body() const {
            std::cout << "\nCelestial body:" << std::endl; 
            print_type(); 
            print_name(); 
            print_mass();
            print_radius();
            print_position();
            std::cout << std::endl; 
        }

};


class Planet : public CelestialBody {

    protected:

        // =============================================
        // class member
        // =============================================
        
        double m_coord_aphelion{}, m_coord_perihelion{}; 
        double m_vel_aphelion{}, m_vel_perihelion{}; 
        double m_period{};

    public: 

        // =============================================
        // constructors and destructor
        // =============================================

        Planet(const char* name) : CelestialBody(name, "Planet") {
            Mass::set_mass_um_prefix("k");
            Position::set_position(zeros(2, 3)); 
            Position::set_coord_um_prefix("k"); 
            Position::set_vel_um_prefix("k"); 
            Sphere::set_center(zeros(3));
            Sphere::set_center_um_prefix("k");

            if (name == "Sun") {
                Mass::set_mass(1.98844E30);
                Sphere::set_radius(695700);
            }

            if (name == "Mercury") {
                Mass::set_mass(0.33010E24);
                Sphere::set_radius(2440.5);
                m_coord_aphelion = 69.818E6;
                m_coord_perihelion = 46E6;
                m_vel_aphelion = 38.86;
                m_vel_perihelion = 58.98;
                m_period = 87.969;
            }
    
            if (name == "Venus") {
                Mass::set_mass(4.8673E24);
                Sphere::set_radius(6051.8); 
                m_coord_aphelion = 108.941E6;
                m_coord_perihelion = 107.480E6;
                m_vel_aphelion = 34.79;
                m_vel_perihelion = 35.26;
                m_period = 224.701;
            }       

            if (name == "Earth") {
                Mass::set_mass(5.9722E24);
                Sphere::set_radius(6378.137);
                m_coord_aphelion = 152.100E6; 
                m_coord_perihelion = 147.095E6; 
                m_vel_aphelion = 29.2911; 
                m_vel_perihelion = 30.2865;
                m_period = 365.256;
            }
            
            if (name == "Mars") {
                Mass::set_mass(0.64169E24);
                Sphere::set_radius(3396.2);
                m_coord_aphelion = 249.261E6;
                m_coord_perihelion = 206.650E6;
                m_vel_aphelion = 21.97;
                m_vel_perihelion = 26.50;
                m_period = 686.980;
            }  

            if (name == "Jupiter") {
                Mass::set_mass(1898.13E24);
                Sphere::set_radius(71492);
                m_coord_aphelion = 816.363E6;
                m_coord_perihelion = 740.595E6;
                m_vel_aphelion = 12.44;
                m_vel_perihelion = 13.72;
                m_period = 4332.589;
            }  

            if (name == "Saturn") {
                Mass::set_mass(568.32E24);
                Sphere::set_radius(60268);
                m_coord_aphelion = 1506.527E6;
                m_coord_perihelion = 1357.554E6;
                m_vel_aphelion = 9.09;
                m_vel_perihelion = 10.18;
                m_period = 10759.22;
            }  

            if (name == "Uranus") {
                Mass::set_mass(86.811E24);
                Sphere::set_radius(25559);
                m_coord_aphelion = 3001.390E6;
                m_coord_perihelion = 2732.696E6;
                m_vel_aphelion = 6.49;
                m_vel_perihelion = 7.11;
                m_period = 30685.4;
            }  

            if (name == "Neptune") {
                Mass::set_mass(102.409E24);
                Sphere::set_radius(24764);
                m_coord_aphelion = 4558.857E6;
                m_coord_perihelion = 4471.050E6;
                m_vel_aphelion = 5.37;
                m_vel_perihelion = 5.50;
                m_period = 60189.; 
            }  

        }


        // =============================================
        // get methods
        // =============================================

        double get_coord_aphelion() const { return m_coord_aphelion; }

        double get_coord_perihelion() const { return m_coord_perihelion; }

        double get_vel_aphelion() const { return m_vel_aphelion; }

        double get_vel_perihelion() const { return m_vel_perihelion; }

        double get_period() const { return m_period; }


        // =============================================
        // print methods
        // =============================================

        void print_period() const { std::cout << "- period = " << get_period() << " days" << std::endl; }

        void print_body() const {
            std::cout << "\nCelestial body:" << std::endl; 
            print_type(); 
            print_name(); 
            print_mass();
            print_radius();
            print_period();
            print_position();
            std::cout << std::endl; 
        }

}; 


// class Satelite : public CelestialBody {

//     protected:

//         // =============================================
//         // class member
//         // =============================================
        
//         double m_coord_apogee{}, m_coord_perigee{}; 
//         double m_vel_apogee{}, m_vel_perigee{}; 
//         const char* m_planet_associated{}; 


//     public:

//         // =============================================
//         // constructors and destructor
//         // =============================================

//         Satelite(const char* name) : CelestialBody(name, "Satelite") {

//             if (name == "Moon") {
//                 set_mass(0.07346e24);
//                 set_radius(1738.1);
//                 m_planet_associated = "Earth";
//                 m_coord_apogee = 0.4055e6;
//                 m_coord_perigee = 0.3633e6;
//                 m_vel_apogee = 0.970;
//                 m_vel_perigee = 1.082;
//             }

//         }
    

//         // =============================================
//         // get methods
//         // =============================================

//         double get_coord_apogee() const { return m_coord_apogee; }

//         double get_coord_perigee() const { return m_coord_perigee; }

//         double get_vel_apogee() const { return m_vel_apogee; }

//         double get_vel_perigee() const { return m_vel_perigee; }

//         const char* get_planet_associated() const { return m_planet_associated; }

            
// }; 



// // class DwarfPlanet : public CelestialBody {

// // };

// // class Luna : public Planet {
// //  public:
// //    Luna() { 
// //       CelestialBody();
// //       set_mass(0.073E24);
// //       set_radius(1737);
// //       setNome("Luna");
// //    } 

// //    double getPosapogee() const { return m_apogee; }
// //    double getPosperigee() const { return m_perigee; }
// //    double getVelapogee() const { return m_vapogee; }
// //    double getVelperigee() const { return m_vperigee; }

// //  private: 
// //    const double m_apogee{0.406E6}, m_perigee{0.363E6}; 
// //    const double m_vapogee{0.971}, m_vperigee{1.083};
// // }; 


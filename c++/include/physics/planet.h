
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Solar System's planets. 
// last updated:    05/07/2022


#include "celestial_body.h"


class Planet : public CelestialBody {

    public:

        // =============================================
        // class member
        // =============================================
        
        double m_coord_afelio{}, m_coord_perielio{}; 
        double m_vel_afelio{}, m_vel_perielio{}; 


        // =============================================
        // constructors and destructor
        // =============================================

        Planet(const char* name) : CelestialBody(name) {

            if (name == "Sun") {
                set_mass(1.98844E30);
                set_radius(695700);
            }

            if (name == "Mercury") { 
                set_mass(0.33010E24);
                set_radius(2440.5);
                m_coord_afelio = 69.818E6;
                m_coord_perielio = 46E6;
                m_vel_afelio = 38.86;
                m_vel_perielio = 58.98;
            }
    
            if (name == "Venus") { 
                set_mass(4.8673E24);
                set_radius(6051.8);
                m_coord_afelio = 108.941E6;
                m_coord_perielio = 107.480E6;
                m_vel_afelio = 34.79;
                m_vel_perielio = 35.26;
            }       

            if (name == "Earth") { 
                set_mass(5.9722E24);
                set_radius(6378.137);
                m_coord_afelio = 152.100E6; 
                m_coord_perielio = 147.095E6; 
                m_vel_afelio = 29.2911; 
                m_vel_perielio = 30.2865;
            }
            
            if (name == "Mars") { 
                set_mass(0.64169E24);
                set_radius(3396.2);
                m_coord_afelio = 249.261E6;
                m_coord_perielio = 206.650E6;
                m_vel_afelio = 21.97;
                m_vel_perielio = 26.50;
            }  

            if (name == "Jupiter") { 
                set_mass(1898.13E24);
                set_radius(71492);
                m_coord_afelio = 816.363E6;
                m_coord_perielio = 740.595E6;
                m_vel_afelio = 12.44;
                m_vel_perielio = 13.72;
            }  

            if (name == "Saturn") { 
                set_mass(568.32E24);
                set_radius(60268);
                m_coord_afelio = 1506.527E6;
                m_coord_perielio = 1357.554E6;
                m_vel_afelio = 9.09;
                m_vel_perielio = 10.18;
            }  

            if (name == "Uranus") { 
                set_mass(86.811E24);
                set_radius(25559);
                m_coord_afelio = 3001.390E6;
                m_coord_perielio = 2732.696E6;
                m_vel_afelio = 6.49;
                m_vel_perielio = 7.11;
            }  

            if (name == "Neptune") { 
                set_mass(102.409E24);
                set_radius(24764);
                m_coord_afelio = 4558.857E6;
                m_coord_perielio = 4471.050E6;
                m_vel_afelio = 5.37;
                m_vel_perielio = 5.50;
            }  

        }

        Planet(const char* name, const std::vector<double>& coord, const std::vector<double>& vel, const double& mass, const double& radius = 0) : 
            CelestialBody(name, coord, vel, mass, radius) {}
        
        Planet(const char* name, const std::vector<std::vector<double>>& pos, const double& mass, const double& radius = 0) : 
            CelestialBody(name, pos, mass, radius) {}


        // =============================================
        // get Afelio & Perielio position methods
        // =============================================

        double get_coord_afelio() const { return m_coord_afelio; }

        double get_coord_perielio() const { return m_coord_perielio; }

        double get_vel_afelio() const { return m_vel_afelio; }

        double get_vel_perielio() const { return m_vel_perielio; }

    
}; 


// class Luna : public Planet {
//  public:
//    Luna() { 
//       CelestialBody();
//       set_mass(0.073E24);
//       set_radius(1737);
//       setNome("Luna");
//    } 

//    double getPosApogeo() const { return m_apogeo; }
//    double getPosPerigeo() const { return m_perigeo; }
//    double getVelApogeo() const { return m_vapogeo; }
//    double getVelPerigeo() const { return m_vperigeo; }

//  private: 
//    const double m_apogeo{0.406E6}, m_perigeo{0.363E6}; 
//    const double m_vapogeo{0.971}, m_vperigeo{1.083};
// }; 


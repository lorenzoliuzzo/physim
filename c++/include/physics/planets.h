
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Basic physics system of generic objects.
// last updated:    03/07/2022


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
            
            if (name == "Earth") { 
                set_mass(5.972E24);
                set_radius(6378);
                m_coord_afelio = 152.098E6; 
                m_coord_perielio = 147.098E6; 
                m_vel_afelio = 29.2911; 
                m_vel_perielio = 30.2865;
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











// class Mercurio : public Planet {
//   public: 
//     Mercurio() { 
//         CelestialBody();
//         setMassa(0.330E24);
//         setRaggio(2439);
//         setNome("Mercurio");
//     }

//     double getAfelio() const { return m_afelio; }
//     double getPerielio() const { return m_perielio; }

//   private:
//     const double m_afelio{69.8E6}, m_perielio{46E6}; 
// }; 

// class Venere : public Planet {
//   public:
//     Venere() { 
//         CelestialBody();
//         setMassa(4.87E24);
//         setRaggio(6502);
//         setNome("Venere");
//     } 

//     double getAfelio() const { return m_afelio; }
//     double getPerielio() const { return m_perielio; }

//   private:
//     const double m_afelio{108.9E6}, m_perielio{107.5E6}; 
// }; 


// class Luna : public Planet {
//  public:
//    Luna() { 
//       CelestialBody();
//       setMassa(0.073E24);
//       setRaggio(1737);
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

// class Marte : public Planet {
//  public:
//    Marte() { 
//       CelestialBody();
//       setMassa(0.642E24);
//       setRaggio(3396);
//       setNome("Marte");
//    } 

//    double getAfelio() const { return m_afelio; }
//    double getPerielio() const { return m_perielio; }

//  private:
//    const double m_afelio{249.3E6}, m_perielio{206.7E6}; 
// }; 

// class Giove : public Planet {
//  public:
//    Giove() { 
//       CelestialBody();
//       setMassa(1898E24);
//       setRaggio(71492);
//       setNome("Giove");
//    } 

//    double getAfelio() const { return m_afelio; }
//    double getPerielio() const { return m_perielio; }

//  private:
//    const double m_afelio{816.4E6}, m_perielio{740.6E6}; 
// }; 

// class Saturno : public Planet {
//  public:
//    Saturno() { 
//       CelestialBody();
//       setMassa(568E24);
//       setRaggio(60268);
//       setNome("Saturno");
//    } 

//    double getAfelio() const { return m_afelio; }
//    double getPerielio() const { return m_perielio; }

//  private:
//    const double m_afelio{1506.5E6}, m_perielio{1357.6E6}; 
// }; 

// class Urano : public Planet {
//  public:
//    Urano() { 
//       CelestialBody();
//       setMassa(86.8E24);
//       setRaggio(25559);
//       setNome("Urano");
//    } 

//    double getAfelio() const { return m_afelio; }
//    double getPerielio() const { return m_perielio; }

//  private:
//    const double m_afelio{3001.4E6}, m_perielio{2732.7E6}; 
// }; 

// class Nettuno : public Planet {
//  public:
//    Nettuno() { 
//       CelestialBody();
//       setMassa(102E24);
//       setRaggio(24764);
//       setNome("Nettuno");
//    } 
//    double getAfelio() const { return m_afelio; }
//    double getPerielio() const { return m_perielio; }

//  private:
//    const double m_afelio{4558.9E6}, m_perielio{4471.1E6}; 
// }; 

// class Plutone : public Planet {
//  public:
//    Plutone() { 
//       CelestialBody();
//       setMassa(0.0130E24);
//       setRaggio(1188);
//       setNome("Plutone");
//    } 

//    double getAfelio() const { return m_afelio; }
//    double getPerielio() const { return m_perielio; }

//  private:
//    const double m_afelio{7375.9E6}, m_perielio{4436.8E6}; 
// }; 

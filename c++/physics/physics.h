
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools. 
// last updated:    13/07/2022


#include <iostream>
#include <cassert>

#include "vector_algebra.h"


namespace physics {

    namespace tools {

        namespace units {

            enum class base_enum { gram = 0, 
                                   meter = 1,
                                   second = 2,
                                   kelvin = 3,
                                   ampere = 4,
                                   mol = 5,  
                                   candela = 6 };    
            
            class base {

                public: 

                    // =============================================
                    // class member
                    // =============================================

                    base_enum m_base; 


                    // =============================================
                    // constructor and destructor
                    // =============================================        
                    
                    base(const int& __base) { m_base = base_enum(__base); }

                    ~base() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    void set_base(const int& __base) { m_base = base_enum(__base); }

                    constexpr base_enum get_base() { return m_base; }


                    // =============================================
                    // print methods
                    // =============================================   
                    
                    void print() const {
                        switch (m_base) {
                            case base_enum::gram:     std::cout << "g" << std::endl;     break;
                            case base_enum::meter:    std::cout << "m" << std::endl;     break;
                            case base_enum::second:   std::cout << "s" << std::endl;     break;
                            case base_enum::kelvin:   std::cout << "K" << std::endl;     break;
                            case base_enum::ampere:   std::cout << "A" << std::endl;     break;
                            case base_enum::mol:      std::cout << "Mol" << std::endl;   break;
                            case base_enum::candela:  std::cout << "cd" << std::endl;    break;
                            default:                  std::cout << "what" << std::endl;  break;
                        }
                    }


                    // =============================================
                    // base methods
                    // =============================================   
                    
                    

            };


            void check_base(base_enum b1, base_enum b2) {
                assert(b1 == b2 && "The bases must be the same! Can't sum bananas and pijamas!");
            }


            enum class prefix_enum { pico = 0,
                                     nano = 1,
                                     micro = 2,
                                     milli = 3,
                                     centi = 4,
                                     deci = 5,
                                     none = 6, 
                                     deca = 7,
                                     hecto = 8,
                                     kilo = 9,
                                     mega = 10,
                                     giga = 11,
                                     tera = 12 };                

            class prefix {

                public: 

                    // =============================================
                    // class member
                    // =============================================

                    prefix_enum m_prefix; 


                    // =============================================
                    // constructor and destructor
                    // =============================================        
                    
                    prefix(const int& __prefix = 6) { m_prefix = prefix_enum(__prefix); }

                    ~prefix() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    void set_prefix(const int& __prefix) { m_prefix = prefix_enum(__prefix); }

                    constexpr prefix_enum get_prefix() { return m_prefix; }


                    // =============================================
                    // print methods
                    // =============================================   
                                                    
                    void print() const {
                        switch (m_prefix) {
                            case prefix_enum::pico:     std::cout << " p";      break;
                            case prefix_enum::nano:     std::cout << " n";      break;
                            case prefix_enum::micro:    std::cout << " micro";  break;
                            case prefix_enum::milli:    std::cout << " m";      break;
                            case prefix_enum::centi:    std::cout << " c";      break;
                            case prefix_enum::deci:     std::cout << " d";      break;
                            case prefix_enum::none:  std::cout << " ";       break;
                            case prefix_enum::deca:     std::cout << " da";     break;
                            case prefix_enum::hecto:    std::cout << " hc";     break;
                            case prefix_enum::kilo:     std::cout << " k";      break;
                            case prefix_enum::mega:     std::cout << " M";      break;
                            case prefix_enum::giga:     std::cout << " G";      break;
                            case prefix_enum::tera:     std::cout << " T";      break;
                            default:                    std::cout << " what";   break;
                        }
                    }        
            
            };


            // void check_prefix(prefix_enum b1, prefix_enum b2) {
            //     if(b1 !== b2) conver(b1, b2);
            // }            

            // void convert(prefix p1, prefix p2) {}


            class unit : public base, prefix {

                public: 

                    // =============================================
                    // constructor and destructor
                    // =============================================

                    unit(const int& __base, const int& __prefix = 6) : base(__base), prefix(__prefix) {}

                    ~unit() {}


                    // =============================================
                    // set and get methods
                    // =============================================
    
                    void set_unit(const int& __base, const int& __prefix) { 
                        base::set_base(__base); 
                        prefix::set_prefix(__prefix);  
                    }

                    unit get_unit() const { return *this; }


                    // =============================================
                    // print methods
                    // =============================================   
                    
                    void print() const {
                        prefix::print(); 
                        base::print(); 
                    }

            };
            
        }

        namespace measurement {
            class measure {

                public:

                    // =============================================
                    // class member
                    // =============================================
                    
                    double m_value;

                    double m_error;  


                    // =============================================
                    // constructor and destructor
                    // =============================================

                    measure(const double& value, const double& error = 0.) : m_value{value}, m_error{error} {}

                    ~measure() {}


                    // =============================================
                    // set and get methods
                    // =============================================       

                    void set_value(const double& value) { m_value = value; }         

                    void set_error(const double& error) { m_error = error; }         

                    constexpr double get_value() const { return m_value; }
                    
                    constexpr double get_error() const { return m_error; }

                    measure get_measure() const { return *this; }

                    
                    // =============================================
                    // print methods
                    // =============================================   

                    void print() const { std::cout << get_value() << " ± " << get_error(); }

            }; 

            class measurement : public measure, public units::unit {

                public:

                    // =============================================
                    // constructor and destructor
                    // =============================================

                    measurement(const double& value, const int& __base, const int& __prefix = 6) : units::unit(__base, __prefix), measure(value) {}

                    measurement(const double& value, const double& error, const int& __base, const int& __prefix = 6) : measure(value, error), units::unit(__base, __prefix) {}

                    ~measurement() {}


                    // =============================================
                    // set and get methods
                    // =============================================       
                    
                    measurement get_measurement() const { return *this; }


                    // // =============================================
                    // // print methods
                    // // =============================================   

                    void print() const {
                        measure::print(); 
                        units::unit::print(); 
                    }

            };

        }


        namespace position {

            class coordinates : public tools::units::unit {

                protected: 

                    // =============================================
                    // class members
                    // =============================================
                
                    // coordinates:     [x] [y] [z] 
                    
                    std::vector<double> m_coordinates = zeros(3);


                public:  

                    // =============================================
                    // constructor and destructor
                    // =============================================

                    coordinates(const std::vector<double>& coord, const int& __prefix = 6) : m_coordinates{coord}, tools::units::unit(1, __prefix) {}

                    ~coordinates() {}
                    
                    // =============================================
                    // set methods
                    // =============================================

                    void set_coordinates(const std::vector<double>& coord) { m_coordinates = coord; }
                    
                    void set_coordinate_x(const double& x) { m_coordinates[0] = x; }

                    void set_coordinate_y(const double& y) { m_coordinates[1] = y;  }

                    void set_coordinate_z(const double& z) { m_coordinates[2] = z; }
                
                    
                    // =============================================
                    // get methods
                    // =============================================

                    std::vector<double> get_coordinates() const { return m_coordinates; }

                    double get_coordinate_x() const { return m_coordinates[0]; }

                    double get_coordinate_y() const { return m_coordinates[1]; }

                    double get_coordinate_z() const { return m_coordinates[2]; }
                    
                    double get_magnitude() const {
                        return sqrt(pow(m_coordinates[0], 2) +                 
                                    pow(m_coordinates[1], 2) + 
                                    pow(m_coordinates[2], 2));
                    }        

                    double get_distance(const std::vector<double>& coord) const {        
                        return sqrt(pow(coord[0] - m_coordinates[0], 2) + 
                                    pow(coord[1] - m_coordinates[1], 2) + 
                                    pow(coord[2] - m_coordinates[2], 2)); 
                    }
                    
                    double get_rho() const { return sqrt(pow(m_coordinates[0], 2) + pow(m_coordinates[1], 2)); }

                    double get_phi() const { return atan2(m_coordinates[1], m_coordinates[0]); }     

                    double get_phi(const std::vector<double>& coord) const { return atan2(coord[1] - m_coordinates[1], coord[0] - m_coordinates[0]); }

                    double get_theta() const { return acos(m_coordinates[2] / get_magnitude()); }
            
                    double get_theta(const std::vector<double>& coord) { return acos((coord[2] - m_coordinates[2]) / get_distance(coord)); }

                    std::vector<double> get_direction() const {
                        return {cos(get_phi()), sin(get_phi()), m_coordinates[2] / get_magnitude()};
                    } 

                    std::vector<double> get_direction(const std::vector<double>& coord1) const {
                        return {cos(get_phi(coord1)), sin(get_phi(coord1)), (coord1[2] - m_coordinates[2]) / get_distance(coord1)};
                    } 
                    

                    // =============================================
                    // print methods
                    // =============================================

                    void print() const {
                        std::cout << "- coordinates = ";
                        for (auto i : m_coordinates) std::cout << "[" << i << "]\t"; 
                        unit::print();
                    }

            };

            class velocity {

                protected: 

                    // =============================================
                    // class members
                    // =============================================
                
                    // velocity:     [x] [y] [z] 
                    
                    std::vector<double> m_velocity = zeros(3);
              

                public:  

                    // =============================================
                    // constructors
                    // =============================================

                    velocity(const std::vector<double>& vel, const int& __prefix = 6) : m_velocity{vel}, tools::units::unit(1, __prefix) {}
                    
                    // =============================================
                    // set methods
                    // =============================================

                    void set_velocity(const std::vector<double>& vel) { m_velocity = vel; }

                    void set_velocity_x(const double& x) { m_velocity[0] = x; }

                    void set_velocity_y(const double& y) { m_velocity[1] = y;  }

                    void set_velocity_z(const double& z) { m_velocity[2] = z; }


                    // =============================================
                    // get methods
                    // =============================================

                    std::vector<double> get_velocity() const { return m_velocity; }

                    double get_velocity_x() const { return m_velocity[0]; }

                    double get_velocity_y() const { return m_velocity[1]; }

                    double get_velocity_z() const { return m_velocity[2]; }
                    
                    double get_magnitude() const {
                        return sqrt(pow(m_velocity[0], 2) +                 
                                    pow(m_velocity[1], 2) + 
                                    pow(m_velocity[2], 2));
                    }        

                    double get_phi() const { return atan2(m_velocity[1], m_velocity[0]); }     

                    double get_theta() const { return acos(m_velocity[2] / get_magnitude()); }
            
                    std::vector<double> get_direction() const {
                        return {cos(get_phi()), sin(get_phi()), cos(get_theta())};
                    } 
                

                    // =============================================
                    // print methods
                    // =============================================

                    void print_velocity() const {
                        std::cout << "- velocity =    ";
                        for (auto i : m_velocity) std::cout << "[" << i << "]\t"; 
                        unit::print();
                    }

            };

        } 

    } 

}

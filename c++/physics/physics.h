
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools. 
// last updated:    13/07/2022


#include <iostream>
#include <cassert>
#include <utility>

#include "vector_algebra.h"


/* ######## NAMESPACE PHYSICS MAP ########

namespace physics
    |
    '----> namespace tools 
            |
            '----> namespace units 
            |       |
            |       '---> enum class base_enum
            |       '---> class base 
            |       '---> enum class prefix_enum
            |       '---> class prefix 
            |
            '----> namespace measurement
            |
            '----> namespace position
        
*/


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

                protected: 
                    
                    const char* int_to_const_chars(int n) const { 
                        if      (n == 0) { return "g";   }
                        else if (n == 1) { return "m";   }
                        else if (n == 2) { return "s";   }
                        else if (n == 3) { return "K";   }
                        else if (n == 4) { return "A";   }
                        else if (n == 5) { return "mol"; }
                        else if (n == 6) { return "cd";  }
                        else {
                            std::cerr << "Invalid int for the base_enum convertion to const char*" << std::endl; 
                            exit(-11); 
                        }
                    }

                    constexpr int const_chars_to_int(const char* __base) const { 
                        if      (__base == "g")   { return 0; } 
                        else if (__base == "m")   { return 1; } 
                        else if (__base == "s")   { return 2; } 
                        else if (__base == "K")   { return 3; } 
                        else if (__base == "A")   { return 4; } 
                        else if (__base == "mol") { return 5; } 
                        else if (__base == "cd")  { return 6; } 
                        else {
                            std::cerr << "Invalid const char* for the base_enum convertion to int" << std::endl; 
                            exit(-11); 
                        }
                    }


                public: 

                    // =============================================
                    // class member
                    // =============================================

                    base_enum m_base; 


                    // =============================================
                    // constructors and destructor
                    // =============================================        
                    
                    base(const char* __base) { m_base = base_enum(const_chars_to_int(__base)); }

                    base(const int& __base) { m_base = base_enum(__base); }

                    ~base() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    void set_base(const int& __base) { m_base = base_enum(__base); }

                    constexpr base_enum get_base() { return m_base; }

                    constexpr int get_base() const { return static_cast<std::underlying_type<base_enum>::type>(m_base); }


                    // =============================================
                    // print methods
                    // =============================================   
                    
                    void print() const { 
                        std::cout << int_to_const_chars(get_base()) << std::endl; 
                    }

            };


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

                private:       
                                    
                    const char* int_to_const_chars(const int& n) const { 
                        if      (n == 0)  { return "p";  }
                        else if (n == 1)  { return "n";  }
                        else if (n == 2)  { return "µ"; }
                        else if (n == 3)  { return "m"; }
                        else if (n == 4)  { return "c"; }
                        else if (n == 5)  { return "d";  }
                        else if (n == 6)  { return "";      }
                        else if (n == 7)  { return "dec";  }
                        else if (n == 8)  { return "hec"; }
                        else if (n == 9)  { return "k";  }
                        else if (n == 10) { return "M";  } 
                        else if (n == 11) { return "G";  } 
                        else if (n == 12) { return "T";  } 
                        else {
                            std::cerr << "Invalid int for the prefix_enum convertion to const char*" << std::endl; 
                            exit(-11); 
                        }
                    }

                    constexpr int const_chars_to_int(const char* __prefix) const { 
                        if      (__prefix == "pico")  { return 0;  }
                        else if (__prefix == "nano")  { return 1;  }
                        else if (__prefix == "micro") { return 2;  }
                        else if (__prefix == "milli") { return 3;  }
                        else if (__prefix == "centi") { return 4;  }
                        else if (__prefix == "deci")  { return 5;  }
                        else if (__prefix == "")      { return 6;  }
                        else if (__prefix == "deca")  { return 7;  }
                        else if (__prefix == "hecto") { return 8;  }
                        else if (__prefix == "kilo")  { return 9;  }
                        else if (__prefix == "mega")  { return 10; }
                        else if (__prefix == "giga")  { return 11; }
                        else if (__prefix == "tera")  { return 12; }
                        else {
                            std::cerr << "Invalid const char* for the prefix_enum convertion to int" << std::endl; 
                            exit(-11);
                        }
                    }


                public: 

                    // =============================================
                    // class member
                    // =============================================

                    prefix_enum m_prefix; 


                    // =============================================
                    // constructors and destructor
                    // =============================================  

                    prefix(const char* __prefix = "") { m_prefix = prefix_enum(const_chars_to_int(__prefix)); }
                                        
                    prefix(const int& __prefix = 6) { m_prefix = prefix_enum(__prefix); }

                    ~prefix() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    void set_prefix(const int& __prefix) { m_prefix = prefix_enum(__prefix); }

                    constexpr prefix_enum get_prefix() { return m_prefix; }

                    constexpr int get_prefix() const { return static_cast<std::underlying_type<prefix_enum>::type>(m_prefix); }


                    // =============================================
                    // print methods
                    // =============================================   
                    
                    void print() const { 
                        std::cout << " " << int_to_const_chars(get_prefix()); 
                    }
     
            
            };


            void check_base(const base_enum& b1, const base_enum& b2) {
                assert(b1 == b2 && "The bases must be the same! Can't sum bananas and pijamas!");
            }

            // void check_prefix(prefix_enum b1, prefix_enum b2) {
            //     if(b1 !== b2) conver(b1, b2);
            // }            

            // void convert(prefix p1, prefix p2) {}

            // void pow(const base_enum& b, const int& n) {
            //     b.
            // }
            

            
            class unit : public base, prefix {

                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================

                    unit(const int& __base, const int& __prefix = 6) : base(__base), prefix(__prefix) {}
                    
                    unit(const char* __base, const char* __prefix = "") : base(__base), prefix(__prefix) {}
                    
                    unit(const int& __base, const char* __prefix = "") : base(__base), prefix(__prefix) {}
                    
                    ~unit() {}


                    // =============================================
                    // set and get methods
                    // =============================================
    
                    void set_unit(const int& __base, const int& __prefix) { 
                        base::set_base(__base); 
                        prefix::set_prefix(__prefix);  
                    }

                    const unit get_unit() const { return *this; }


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
                    // constructors and destructor
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
                    // constructors and destructor
                    // =============================================

                    measurement(const double& value, const int& __base, const int& __prefix = 6) : units::unit(__base, __prefix), measure(value) {}
                    
                    measurement(const double& value, const char* __base, const char* __prefix = "") : units::unit(__base, __prefix), measure(value) {}

                    measurement(const double& value, const double& error, const int& __base, const int& __prefix = 6) : measure(value, error), units::unit(__base, __prefix) {}
                    
                    measurement(const double& value, const double& error, const char* __base, const char* __prefix = "") : measure(value, error), units::unit(__base, __prefix) {}

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
                    // constructors and destructor
                    // =============================================

                    coordinates(const std::vector<double>& coord, const int& __prefix = 6) : m_coordinates{coord}, tools::units::unit(1, __prefix) {}
                    
                    coordinates(const std::vector<double>& coord, const char* __prefix = "") : m_coordinates{coord}, tools::units::unit(1, __prefix) {}

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

            // class velocity {

            //     protected: 

            //         // =============================================
            //         // class members
            //         // =============================================
                
            //         // velocity:     [x] [y] [z] 
                    
            //         std::vector<double> m_velocity = zeros(3);
              

            //     public:  

            //         // =============================================
            //         // constructors
            //         // =============================================

            //         velocity(const std::vector<double>& vel, const int& __prefix = 6) : m_velocity{vel}, tools::units::unit(1, __prefix) {}
                    
            //         // =============================================
            //         // set methods
            //         // =============================================

            //         void set_velocity(const std::vector<double>& vel) { m_velocity = vel; }

            //         void set_velocity_x(const double& x) { m_velocity[0] = x; }

            //         void set_velocity_y(const double& y) { m_velocity[1] = y;  }

            //         void set_velocity_z(const double& z) { m_velocity[2] = z; }


            //         // =============================================
            //         // get methods
            //         // =============================================

            //         std::vector<double> get_velocity() const { return m_velocity; }

            //         double get_velocity_x() const { return m_velocity[0]; }

            //         double get_velocity_y() const { return m_velocity[1]; }

            //         double get_velocity_z() const { return m_velocity[2]; }
                    
            //         double get_magnitude() const {
            //             return sqrt(pow(m_velocity[0], 2) +                 
            //                         pow(m_velocity[1], 2) + 
            //                         pow(m_velocity[2], 2));
            //         }        

            //         double get_phi() const { return atan2(m_velocity[1], m_velocity[0]); }     

            //         double get_theta() const { return acos(m_velocity[2] / get_magnitude()); }
            
            //         std::vector<double> get_direction() const {
            //             return {cos(get_phi()), sin(get_phi()), cos(get_theta())};
            //         } 
                

            //         // =============================================
            //         // print methods
            //         // =============================================

            //         void print_velocity() const {
            //             std::cout << "- velocity =    ";
            //             for (auto i : m_velocity) std::cout << "[" << i << "]\t"; 
            //             unit::print();
            //         }

            // };

        } 

    } 

}


// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Mathematical functions
// last updated:    02/07/2022


#pragma once
#include "vector_algebra.h"
#include <iomanip>


class FunctionBase {

    public: 
    
        // =============================================
        // destructor
        // =============================================     
    
        virtual ~FunctionBase() {}

    
        // =============================================
        // eval methods
        // =============================================
    
        virtual double eval(const double& x) const = 0; 
    
    
        // =============================================
        // print methods
        // =============================================
    
        void print_eval(const double& x, const double& precision = 1.e-5) const {
            std::cout << "f (" << x << ") = " << std::fixed << std::setprecision((int)-log10(precision)) << eval(x) << std::endl; 
        }
    
        virtual void print_equation() const = 0;
    
        
        // =============================================
        // extra methods
        // =============================================
    
        int signum(const double& x) const { return (x == 0. ? 0. : (x > 0 ? 1. : -1)); }
    
}; 


class Square : public FunctionBase {

    private: 

        // =============================================
        // class members
        // =============================================
        
        double m_a, m_b, m_c, m_delta; 


    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================     
  
        Square(double a, double b, double c) : m_a{a}, m_b{b}, m_c{c}, m_delta{pow(m_b, 2) - 4 * m_a * m_c} {}
        
        ~Square() {} 

  
        // =============================================
        // set methods
        // =============================================
  
        void set_a(double a) { m_a = a; } 

        void set_b(double b) { m_b = b; }

        void set_c(double c) { m_c = c; }  
  
  
        // =============================================
        // get methods
        // =============================================

        double get_a() const { return m_a; }

        double get_b() const { return m_b; } 

        double get_c() const { return m_c; } 

        double get_delta() const { return m_delta; } 

        std::vector<double> get_roots() const {
            std::vector<double> appo = zeros(2); 
            if (m_delta == 0) { 
                appo[0] = - m_b / (2 * m_a); 
                appo[1] = - m_b / (2 * m_a); 
            } else if (m_delta > 0) {
                appo[0] = (- m_b - sqrt(m_delta)) / (2 * m_a); 
                appo[1] = (- m_b + sqrt(m_delta)) / (2 * m_a);
            } else std::cout << std::endl << "There are not real solutions..." << std::endl << std::endl; // note: create a complex class then come back here
            return appo; 
        }


        // =============================================
        // eval methods
        // =============================================

        double eval_Horner(double x) const { return m_c + x * (m_b + x * m_a); }

        double eval(const double& x) const override { return m_a * pow(x, 2) + m_b * x + m_c; }
  
  
        // =============================================
        // print methods
        // =============================================

        void print_equation() const override {
            std::cout << "Equation:  y = " << m_a << "x^2 + " << m_b << "x + " << m_c << std::endl;
        }
    
        void print_roots(double precision = 1.e-5) const {
            std::vector<double> appo = get_roots(); 
            std::cout << "Roots: x = " << std::fixed << std::setprecision((int)-log10(precision)) << appo[0] << "  V  x = " << appo[1] << std::endl; 
        }

}; 


class Cube : public FunctionBase {
    
    private: 

        // =============================================
        // class members
        // =============================================
        
        double m_a, m_b, m_c, m_d; 


    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
    
        Cube(double a, double b, double c, double d) : m_a{a}, m_b{b}, m_c{c}, m_d{d} {}
    
        ~Cube() {} 
    
    
        // =============================================
        // set methods
        // =============================================
  
        void set_a(double a) { m_a = a; } 

        void set_b(double b) { m_b = b; }

        void set_c(double c) { m_c = c; }  
  
        void set_d(double d) { m_d = d; } 
   
    
        // =============================================
        // get methods
        // =============================================

        double get_a() const { return m_a; }

        double get_b() const { return m_b; } 

        double get_c() const { return m_c; } 

        double get_d() const { return m_d; } 

    
        // =============================================
        // eval methods
        // =============================================

        double eval(const double& x) const override { return m_a * pow(x, 3) + m_b * pow(x, 2) + m_c * x + m_d; }
  
  
        // =============================================
        // print methods
        // =============================================

        void print_equation() const override {
            std::cout << "Equation:  y = " << m_a << "x ^ 3 + " << m_b << "x ^ 2 + " << m_c << "x + " << m_d << std::endl;
        }
  
}; 


class SquareRoot : public FunctionBase {
    
    private: 

        // =============================================
        // class members
        // =============================================
        
        // y = c1 * x ^ (1 / 2)
        double m_c1; 


    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
    
        SquareRoot(double c1 = 1) : m_c1{c1} {}
    
        ~SquareRoot() {} 
    
    
        // =============================================
        // set & get methods
        // =============================================
    
        void set_c1(double c1) { m_c1 = c1; }

        double get_c1() const { return m_c1; }

    
        // =============================================
        // eval methods
        // =============================================

        double eval(const double& x) const override { return m_c1 * pow(x, 0.5); }
    
        double eval(const FunctionBase& f, double x) const { return m_c1 * pow(f.eval(x), 0.5); }
  
  
        // =============================================
        // print methods
        // =============================================

        void print_equation() const override {
            std::cout << "Equation:  y = " << m_c1 << "x ^ (1/2)" << std::endl;
        }        
    
    
};


class CubeRoot : public FunctionBase {
    
    private: 

        // =============================================
        // class members
        // =============================================
        
        // y = c1 * x ^ (1 / 3)
        double m_c1; 


    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
    
        CubeRoot(double c1 = 1) : m_c1{c1} {}
    
        ~CubeRoot() {} 
    
    
        // =============================================
        // set & get methods
        // =============================================
    
        void set_c1(double c1) { m_c1 = c1; }

        double get_c1() const { return m_c1; }

    
        // =============================================
        // eval methods
        // =============================================

        double eval(const double& x) const override { return m_c1 * pow(x, 1. / 3.); }
    
        double eval(const FunctionBase& f, double x) const { return m_c1 * pow(f.eval(x), 1. / 3.); }
  
  
        // =============================================
        // print methods
        // =============================================

        void print_equation() const override {
            std::cout << "Equation:  y = " << m_c1 << "x^(1/3)" << std::endl;
        }        
    
};


class Exponential : public FunctionBase {
  
    private:
  
        // =============================================
        // class members
        // =============================================
        
        // y = c1 * base ^ (c2 * x)
        double m_base; 
        double m_c1, m_c2;

    
    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
  
        Exponential(double base, double c1 = 1, double c2 = 1) : m_base{base}, m_c1{c1}, m_c2{c2} {}
    
        ~Exponential() {} 
  

        // =============================================
        // set methods
        // =============================================
  
        void set_base(double base) { m_base = base; } 
    
        void set_c1(double c1) { m_c1 = c1; }
    
        void set_c2(double c2) { m_c2 = c2; }
      
  
        // =============================================
        // get methods
        // =============================================
  
        double get_base() { return m_base; }
  
        double get_c1() { return m_c1; }

        double get_c2() { return m_c2; }
    
    
        // =============================================
        // eval methods
        // =============================================
  
        double eval(const double& x) const override { return m_c1 * pow(m_base, m_c2 * x); }
  
        double eval(const FunctionBase& f, double x) const { return m_c1 * pow(m_base, m_c2 * f.eval(x)); }
    
    
        // =============================================
        // print methods
        // =============================================
  
        void print_equation() const override {
            std::cout << "Equation:  y = " << m_c1 << " * " << m_base << "^ (" << m_c2 << " * x) " << std::endl;
        }

}; 


class Logarithm : public FunctionBase {
  
    private:
  
        // =============================================
        // class members
        // =============================================
        
        // y = c1 * log_base (c2 * x) 
        double m_base; 
        double m_c1, m_c2; 

    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
  
        Logarithm(double base, double c1 = 1, double c2 = 1) : m_base{base}, m_c1{c1}, m_c2{c2} {}
    
        ~Logarithm() {} 
  

        // =============================================
        // set methods
        // =============================================
  
        void set_base(double base) { m_base = base; } 
  
        void set_c1(double c1) { m_c1 = c1; }
    
        void set_c2(double c2) { m_c2 = c2; }
    
    
        // =============================================
        // get methods
        // =============================================
  
        double get_base() { return m_base; }

        double get_c1() { return m_c1; }

        double get_c2() { return m_c2; }
    

        // =============================================
        // eval methods
        // =============================================
  
        double eval(const double& x) const override { return m_c1 * log(m_c2 * x) / log(m_base); }
        
        double eval(const FunctionBase& f, double x) const { return m_c1 * log(m_c2 * f.eval(x)) / log(m_base); }
  
    
        // =============================================
        // print methods
        // =============================================
  
        void print_equation() const override {
            std::cout << "Equation:  y = " << m_c1 << " * log_ " << m_base << "(" << m_c2 << " * x)" << std::endl;
        }

}; 


class Sine : public FunctionBase {
    
    private: 
        // =============================================
        // class members
        // =============================================   
        
        // y = c1 * sin (c2 * x)
        double m_c1, m_c2; 
    
    
    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
  
        Sine(double c1 = 1, double c2 = 1) : m_c1{c1}, m_c2{c2} {} 
    
        ~Sine() {} 

    
        // =============================================
        // set methods
        // =============================================
  
        void set_c1(double c1) { m_c1 = c1; }
    
        void set_c2(double c2) { m_c2 = c2; }
    
    
        // =============================================
        // get methods
        // =============================================

        double get_c1() { return m_c1; }

        double get_c2() { return m_c2; }
    
    
        // =============================================
        // eval methods
        // =============================================
  
        double eval(const double& x) const override { return m_c1 * sin(m_c2 * x); }
  
        double eval(const FunctionBase& f, double x) const { return m_c1 * sin(m_c2 * f.eval(x)); }
  
    
        // =============================================
        // print methods
        // =============================================
  
        void print_equation() const override {
            std::cout << "Equation:  y = sin(x)" << std::endl;
        }

}; 


class Cosine : public FunctionBase {
    
    private: 
        // =============================================
        // class members
        // =============================================   
        
        // y = c1 * cos (c2 * x)
        double m_c1, m_c2; 
    
    
    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
  
  
        Cosine(double c1 = 1, double c2 = 1) : m_c1{c1}, m_c2{c2} {} 
    
        ~Cosine() {} 

    
        // =============================================
        // set methods
        // =============================================
  
        void set_c1(double c1) { m_c1 = c1; }
    
        void set_c2(double c2) { m_c2 = c2; }
    
    
        // =============================================
        // get methods
        // =============================================

        double get_c1() { return m_c1; }

        double get_c2() { return m_c2; }
    

        // =============================================
        // eval methods
        // =============================================
  
        double eval(const double& x) const override { return m_c1 * cos(m_c2 * x); }
  
        double eval(const FunctionBase& f, double x) const { return m_c1 * cos(m_c2 * f.eval(x)); }
  
  
        // =============================================
        // print methods
        // =============================================
  
        void print_equation() const override {
            std::cout << "Equation:  y = cos(x)" << std::endl;
        }

}; 


// class FunctionSum {

//     private: 
//         std::vector<FunctionBase> m_sum;
//         double m_constant; 

//     public:

//         FunctionSum(std::vector<FunctionBase> functions, const double& c = 0) : m_sum{functions}, m_constant{c} {}

//         ~FunctionSum() {}

//         double eval(const double& x) const {
//             double sum{m_constant}; 
//             for (unsigned int i{}; i < m_sum.size(); i++) sum += m_sum[i].eval(x); 
//             return sum; 
//         }

// }; 


// class FunctionSub {

//     private: 
//         std::vector<FunctionBase> m_sub; 
//         double m_constant; 

//     public:

//         FunctionSub(const std::vector<FunctionBase>& functions, const double& c = 0) : m_sub{functions}, m_constant{c} {}

//         ~FunctionSub() {}

//         double eval(const double& x) {
//             double sub{m_constant}; 
//             for (unsigned int i{}; i < m_sub.size(); i++) sub -= m_sub[i].eval(x); 
//             return sub; 
//         }

// }; 



// class FunctionDot {

//     private: 
//         std::vector<FunctionBase> m_dot;
//         double m_constant; 

//     public:

//         FunctionDot(const std::vector<FunctionBase>& functions, const double& c = 1) : m_dot{functions}, m_constant{c} {}

//         ~FunctionDot() {}

//         double eval(const double& x) {
//             double dot{}; 
//             for (unsigned int i{}; i < m_dot.size(); i++) dot *= m_dot[i].eval(x); 
//             return m_constant * dot; 
//         }

// }; 


// class FunctionDiv {

//     private: 
//         std::vector<FunctionBase> m_div;
//         double m_constant; 

//     public:

//         FunctionDiv(const std::vector<FunctionBase>& functions, const double& c = 1) : m_div{functions}, m_constant{c} {}

//         ~FunctionDiv() {}

//         double eval(const double& x) {
//             double div{m_constant}; // + m_div[0].eval(x)}; 
//             for (unsigned int i{1}; i < m_div.size(); i++) div /= m_div[i].eval(x); 
//             return m_constant * div; 
//         }

// }; 


// class FunctionPow {

//     private: 
//         std::vector<FunctionBase> m_div, 
//         double m_constant; 

//     public:

//         FunctionDiv(const std::vector<FunctionBase>& functions, const double& c = 1) : m_div{functions}, m_constant{c} {}

//         ~FunctionDiv() {}

//         double eval(const double& x) {
//             double div{m_constant + m_div[0].eval(x)}; 
//             for (unsigned int i{1}; i < m_div.size(); i++) div /= m_div[i].eval(x); 
//             return m_constant * div; 
//         }

// }; 


// class FunctionCompose {

//     private: 
//         std::vector<FunctionBase> m_div, 
//         double m_constant; 

//     public:

//         FunctionDiv(const std::vector<FunctionBase>& functions, const double& c = 1) : m_div{functions}, m_constant{c} {}

//         ~FunctionDiv() {}

//         double eval(const double& x) {
//             double div{m_constant + m_div[0].eval(x)}; 
//             for (unsigned int i{1}; i < m_div.size(); i++) div /= m_div[i].eval(x); 
//             return m_constant * div; 
//         }

// }; 


// class Function : public FunctionBase {

//     private: 
//         // =============================================
//         // class members
//         // =============================================   
        
//         std::vector<const FunctionBase&> m_functions; 
    
    
//     public:
         
//         // =============================================
//         // constructor and destructor
//         // =============================================   
  
//         Function(std::vector<const FunctionBase&> f) : m_functions{f} {}
    
//         ~Function() {} 

// };














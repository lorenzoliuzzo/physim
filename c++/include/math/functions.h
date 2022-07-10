
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Mathematical functions
// last updated:    03/07/2022


#pragma once
#include "vector_algebra.h"
#include <iomanip>


class FunctionBase {

    public: 
    

        // =============================================
        // constructor and destructor
        // =============================================     
    
        FunctionBase() {}

        ~FunctionBase() {}

    
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


class Functor : public FunctionBase {

    private:

        // =============================================
        // class members
        // =============================================
    
        char m_op; 
        FunctionBase * m_f; 
        FunctionBase * m_g;
    
    
    public:

        // =============================================
        // constructor and destructor
        // =============================================     
    
        Functor(const char& op, FunctionBase& f, FunctionBase& g) {
            m_f = &f; 
            m_g = &g;
            m_op = op;
        }
        
        ~Functor() {}
    
    
        // =============================================
        // eval methods
        // =============================================
    
        double eval(const double& x) const override {
            switch (m_op) {
                case '+':
                    return m_f->eval(x) + m_g->eval(x);

                case '-': 
                    return m_f->eval(x) - m_g->eval(x);

                case '*':
                    return m_f->eval(x) * m_g->eval(x);

                case '/': 
                    return m_f->eval(x) / m_g->eval(x);
                
                case '^':
                    return pow(m_f->eval(x), m_g->eval(x));

                case 'c':
                    return m_f->eval(m_g->eval(x));
            }
            return NAN;
        }

        
        // =============================================
        // print methods
        // =============================================
    
        void print_equation() const override {}
        
};


class Quadratic : public FunctionBase {

    private: 

        // =============================================
        // class members
        // =============================================
        
        double m_a, m_b, m_c, m_delta; 


    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================     
  
        Quadratic(double a, double b, double c) : m_a{a}, m_b{b}, m_c{c}, m_delta{pow(m_b, 2) - 4 * m_a * m_c} {}
        
        ~Quadratic() {} 

  
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
            exit(-11);  
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


class Cubic : public FunctionBase {
    
    private: 

        // =============================================
        // class members
        // =============================================
        
        double m_a, m_b, m_c, m_d; 


    public: 
  
        // =============================================
        // constructor and destructor
        // =============================================   
    
        Cubic(double a, double b, double c, double d) : m_a{a}, m_b{b}, m_c{c}, m_d{d} {}
    
        ~Cubic() {} 
    
    
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


class CubicRoot : public FunctionBase {
    
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
    
        CubicRoot(double c1 = 1) : m_c1{c1} {}
    
        ~CubicRoot() {} 
    
    
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


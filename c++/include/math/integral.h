
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Integral(class) and its derived classes present are some usefull tools 
//                  for numerically evaluate the integral of a function in a range [a, b].
// last updated:    03/07/2022


#pragma once
#include "functions.h"
#include "statistics.h"
#include "random_generator.h"


class Integral {

    protected: 

        // =============================================
        // class members
        // =============================================

        double m_sum, m_integral, m_old_integral, m_error;  
        double m_a, m_b, m_h;
        unsigned int m_steps;
        int m_sign; 
        RandomGenerator m_rand_gen;


    private: 

        // =============================================
        // set methods
        // =============================================

        void set_a(const double& a) { m_a = a; } 

        void set_b(const double& b) { m_b = b; } 

        void set_range(const double& a, const double& b) { m_a = a; m_b = b; }

        void set_steps(unsigned int n) { m_steps = n; }        
        
        void reset_integral() { m_integral = 0; }    

        void set_h() { m_h = fabs(m_b - m_a) / m_steps; }

        void set_sum(const double& x = 0) { m_sum = x; }


    public: 

        // =============================================
        // constructors
        // =============================================

        Integral() : m_a{}, m_b{}, m_h{}, m_sign{}, m_sum{}, m_integral{}, m_old_integral{}, m_error{} {}
        ~Integral() {}

    
        // =============================================
        // get methods
        // =============================================

        double get_a() const { return m_a; }

        double get_b() const { return m_b; }

        unsigned int get_steps() const { return m_steps; }

        int get_sign() const { return m_sign; }

        double get_h() const { return m_h; }

        double get_sum() const { return m_sum; }

        double get_integral() const { return m_integral; }

        double get_error_integral() const { return m_error; }


        // =============================================
        // usefull methods
        // =============================================

        void check_range(const double& a, const double& b) {
            if (a > b) { 
                set_range(b, a);
                m_sign = -1; 
            } else { 
                set_range(a, b); 
                m_sign = 1; 
            }
        }

        void begin_integration(const double& a, const double& b, unsigned int n = 1000, const double& sum0 = 0) {
            check_range(a, b); 
            reset_integral(); 
            set_steps(n); 
            set_h(); 
            set_sum(sum0); 
        }

        void error_integral(const double& new_old) { m_error = 4 * new_old / 3.; } 

        void print_integral(const double& precision = 1.e-5) {
            std::cout << "\nIntegral of f(x) in [" << m_a << ", " << m_b << "] = " << std::setprecision((int)-log10(precision)) << m_integral << std::endl;
        }

        void print_error() {
            std::cout << "error = " << m_error << std::endl;
        }        
        
        
        // =============================================
        // integration methods
        // =============================================

        void midpoint(const double& a, const double& b, const FunctionBase& f, unsigned int n = 1000) {
            begin_integration(a, b, n); 
            for (unsigned int i{}; i < m_steps; i++) { m_sum += (f.eval(m_a + (i + 0.5) * m_h)); }
            m_integral = m_sign * m_sum * m_h; 
        }

        void midpoint_fixed(const double& a, const double& b, const FunctionBase& f, const double& prec = 1.e-5) {
            begin_integration(a, b, 1); 
            double oldintegral;  
            while (true) {
                oldintegral = m_integral; 
                midpoint(m_a, m_b, f, m_steps * 2);
                error_integral(fabs(m_integral - oldintegral));
                if (m_error < prec) break;
            }    
        }
        
        void trapexoid(const double& a, const double& b, const FunctionBase& f, unsigned int n = 1000) {
            begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 2.);
            for (unsigned int i{1}; i < m_steps; i++) { m_sum += f.eval(m_a + i * m_h); }
            m_integral = m_sign * m_sum * m_h; 
        } 

        void trapexoid_fixed(const double& a, const double& b, const FunctionBase& f, const double& prec = 1.e-5) {
            begin_integration(a, b, 2, f.eval(a) + f.eval(b) / 2. + f.eval((a + b) / 2.)); 
            double oldintegral;  
            while (true) {
                oldintegral = m_integral; 
                trapexoid(m_a, m_b, f, m_steps * 2);
                error_integral(fabs(m_integral - oldintegral));
                if (m_error < prec) break;
            }
        }

        void simpson(const double& a, const double& b, const FunctionBase& f, unsigned int n = 1000) {
            if (n % 2 == 0) begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 3.);
            else begin_integration(a, b, n + 1);  
            for (unsigned int i{1}; i < m_steps; i++) {
                m_sum += 2 * (1 + i % 2) * (f.eval(m_a + i * m_h)) / 3.; 
            }
            m_integral = m_sign * m_sum * m_h; 
        }

        void simpson_fixed(const double& a, const double& b, const FunctionBase& f, const double& prec = 1.e-5) {
            begin_integration(a, b, 2, (f.eval(a) + f.eval(b)) / 3.); 
            double oldintegral; 
            while (true) {
                oldintegral = m_integral; 
                simpson(m_a, m_b, f, m_steps * 2);
                error_integral(fabs(m_integral - oldintegral));
                if (m_error < prec) break; 
            }
        }

        void mean(const double& a, const double& b, const FunctionBase& f, unsigned int n = 1000) {
            begin_integration(a, b, n); 
            for (unsigned int i{}; i < n; i ++) {
                m_sum += f.eval(m_rand_gen.rand(a, b)); 
            }
            m_integral = (b - a) * m_sum / n; 
        }

        void mean_fixed(const double& a, const double& b, const FunctionBase& f, const double& prec = 1.e-5) {
            std::vector<double> k{};
            for (unsigned i{}; i < 10000; i++) {
                mean(a, b, f);
                k.push_back(m_integral); 
            }
            double k_mean = sqrt(100) * sd<double>(k); 
            unsigned int N = (unsigned int) pow(k_mean / prec, 2); 
            mean(a, b, f, N); 
        }
  
        void hit_miss(const double& a, const double& b, const FunctionBase& f, const double& fmax, unsigned int n = 1000) {
            begin_integration(a, b, n); 
            double x{}, y{}; 
            unsigned int hits{};
            for (unsigned int i{}; i < n; i ++) {
                x = m_rand_gen.rand(a, b); 
                y = m_rand_gen.rand(0., fmax);  
                if (y <= f.eval(x)) hits++; 
            }
            m_integral = (b - a) * fmax * hits / n; 
        }

        void hit_miss_fixed(const double& a, const double& b, const FunctionBase& f, const double& fmax, const double& prec = 1.e-5) {
            std::vector<double> k{};
            for (unsigned i{}; i < 10000; i++) {
                hit_miss(a, b, f, fmax);
                k.push_back(m_integral); 
            }
            double k_mean = sqrt(100) * sd<double>(k); 
            unsigned int N = (unsigned int) pow(k_mean / prec, 2); 
            hit_miss(a, b, f, fmax, N); 
        }


}; 



// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     RandomGenerator(class) can generate pseudo-casual numbers, following some of the most common distributions of probability.
// last updated:    02/07/2022


#pragma once
#include <iostream>
#include <cmath>


class RandomGenerator {
    
    private: 
    
        // =============================================
        // class members
        // =============================================
    
        unsigned int m_a, m_c, m_m, m_seed;

    
    public: 
  
        // =============================================
        // constructors
        // =============================================
  
        RandomGenerator() { set_a(1664525); set_c(1013904223); set_m(pow(2, 31)); }
  
        RandomGenerator(unsigned int seed) : RandomGenerator() { m_seed = seed; }

  
        // =============================================
        // set methods
        // =============================================

        void set_a(unsigned int a) { m_a = a; }
  
        void set_c(unsigned int c) { m_c = c; }
  
        void set_m(unsigned int m) { m_m = m; }
  
        void set_seed(unsigned int seed) { m_seed = seed; }

  
        // =============================================
        // get methods
        // =============================================

        unsigned int get_a() const { return m_a; }
  
        unsigned int get_c() const { return m_c; }
  
        unsigned int get_m() const { return m_m; }
  
        unsigned int get_seed() const { return m_seed; }
  
  
        // =============================================
        // distributions methods
        // =============================================
  
        double rand(double min = 0., double max = 1.) {
            set_seed((unsigned int)((get_a() * get_seed() + get_c()) % get_m())); 
            return min + (max - min) * get_seed() / get_m(); 
        }

        double exp(double mean) {
            return - log(1 - rand()) / mean; 
        }

        double gauss_box_muller(double mean, double sigma) {
            double s{rand()}, t{rand()}, x{};
            x = sqrt(-2 * log(s)) * cos(2 * M_PI * t);
            return mean + sigma * x;
        }

        double gauss_accept_reject(double mean, double sigma) {
            double x{}, y{}, g{}; 
            while (true) {
                x = rand(-5., 5.); 
                y = rand(); 
                g = exp(- pow(x, 2) / 2); 
                if (y <= g) break;
            }
            std::cout << "sium\n";
            return mean + x * sigma;
        }

}; 

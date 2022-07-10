
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Simple Time(class) and Timer(class) for keeping track of the inexorable passing of time.
// last updated:    10/07/2022


#pragma once
#include <iostream>
#include <chrono>


class Time {

    public:     

        // =============================================
        // class members
        // =============================================     

        double m_time; 

        const char* m_time_um{"s"}; 
        
        const char* m_time_um_prefix{}; 


        // =============================================
        // set and get methods
        // =============================================
        
        void set_time_um_prefix(const char* um_prefix) { m_time_um_prefix = um_prefix; }

        const char* get_time_um() const { return m_time_um; }

        const char* get_time_um_prefix() const { return m_time_um_prefix; }


        // =============================================
        // constructor and destructor
        // =============================================   

        Time(double time = 0, const char* um_prefix = "") : m_time{time}, m_time_um_prefix{um_prefix} {}

        ~Time() {}

        
        // =============================================
        // time methods
        // =============================================   

        double get_time() const { return m_time; }

        void increase_time(const double& h) { m_time += h; } 

        void reset_time() { m_time = 0; }

}; 


class Timer : public Time {

    public:

        // =============================================
        // class members
        // =============================================     
        
        std::chrono::duration<double> m_elapsed_seconds;
        std::chrono::time_point<std::chrono::system_clock> m_start, m_pause, m_end;
        std::time_t m_end_time;


        // =============================================
        // constructor and destructor
        // =============================================   

        Timer(const char* um_prefix = "") : Time(0., um_prefix) {}

        ~Timer() {}

        
        // =============================================
        // timer methods
        // =============================================   

        void start() { m_start = std::chrono::system_clock::now(); }

        void pause() { m_pause = std::chrono::system_clock::now(); }

        void end() { m_end = std::chrono::system_clock::now(); }

        double get_elapsed_time() { 
            m_elapsed_seconds = m_pause - m_start;
            return m_elapsed_seconds.count(); 
        }

        double get_duration() { 
            m_elapsed_seconds = m_end - m_start;
            m_end_time = std::chrono::system_clock::to_time_t(m_end);
            return m_elapsed_seconds.count(); 
        }

        
        // =============================================
        // print time methods
        // =============================================   
        
        void print_elapsed_time() { std::cout << "Elapsed time = " << get_elapsed_time() << " " << get_time_um_prefix() << get_time_um() << std::endl; }

        void print_duration() { std::cout << "Duration time = " << get_duration() << " " << get_time_um_prefix() << get_time_um() << std::endl; }

}; 


// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Basic physics system of generic objects.
// last updated:    05/07/2022


#include "../math/vector_algebra.h"


template <typename T> 
class System {

    private:
        
        // =============================================
        // class members
        // =============================================     

        std::vector<T> m_objects; 
        unsigned int m_objects_count; 


    public:

        // =============================================
        // constructor and destructor
        // =============================================     

        System() {}

        ~System() {}


        // =============================================
        // set & get methods
        // =============================================     
        
        void add_object(T x) { m_objects.push_back(x); } 

        void reset_objects() {
            m_objects.clear(); 
        }

        int get_objects_count() const {
            return m_objects.size();
        }

        std::vector<T> &get_objects() {
            return m_objects;
        }

        T &get_object(unsigned pos) {
            return m_objects[pos]; 
        }

};

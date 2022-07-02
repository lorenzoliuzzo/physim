

struct box {
    
    double m_border_inf_x, m_border_inf_y, m_border_inf_z; 
    double m_border_sup_x, m_border_sup_y, m_border_sup_z; 
    bool m_collision; 

}; 

bool check_collision(const Position& p1, const Position& )

    private: 

        // =============================================
        // class member
        // =============================================

        double m_border_inf_x, m_border_inf_y, m_border_inf_z; 
        double m_border_sup_x, m_border_sup_y, m_border_sup_z; 
        bool m_collision; 


        // =============================================
        // set member
        // =============================================

        void set_border_x(double x_inf, double x_sup) {
            m_bord_inf_x = x_inf; 
            m_bord_sup_x = x_sup; 
        }

        void set_border_y(double y_inf, double y_sup) {
            m_bord_inf_y = y_inf; 
            m_bord_sup_y = y_sup; 
        }

        void set_border_z(double z_inf, double z_sup) {
            m_bord_inf_z = z_inf; 
            m_bord_sup_z = z_sup; 
        }



}; 



void set_precision(double prec) { m_precision = prec; }

double angle_scattering() {}

bool check_collision() { return m_distance}


// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Shape(class) mother class of the most common shapes in the universe.
// last updated:    10/07/2022


#pragma once
#include "coordinates.h"


class Shape3D {

    public: 

        // =============================================
        // class member
        // =============================================

        Coordinates m_center; 

        double m_area, m_volume; 

        const char* m_shape_um{"m"}; 

        const char* m_shape_um_prefix; 


        // =============================================
        // virtual destructor
        // =============================================

        virtual ~Shape3D() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_center(const std::vector<double>& coord) { m_center.set_coordinates(coord); }

        void set_center_um_prefix(const char* um_prefix) { m_shape_um_prefix = um_prefix; }
        
        std::vector<double> get_center() const { return m_center.get_coordinates(); }

        double get_area() const { return m_area; }

        double get_volume() const { return m_volume; }

        const char* get_shape_um() const { return m_shape_um; }
        
        const char* get_shape_um_prefix() const { return m_shape_um_prefix; }
        

        // =============================================
        // print methods
        // =============================================

        void print_area() const { std::cout << "- area = " << get_area() << " (" << get_shape_um_prefix() << get_shape_um() << ") ^ 2 " << std::endl; }

        void print_volume() const { std::cout << "- volume = " << get_volume() << " (" << get_shape_um_prefix() << get_shape_um() << ") ^ 3 " << std::endl; }

}; 


class Cylinder : public Shape3D {

	public:
        
        // =============================================
        // class member
        // =============================================

        double m_radius; 

        double m_height;


        // =============================================
        // constructor and destructor
        // =============================================

		Cylinder(const double& radius, const double& height, const std::vector<double> center, const char* shape_um_prefix = "") : 
            m_radius{radius}, m_height{height} {
            m_center = center;
            m_area = 2 * M_PI * (radius * height +  pow(radius, 2));
            m_volume = 2 * M_PI * radius * height;
            m_shape_um_prefix = shape_um_prefix;  
        }

        ~Cylinder() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_radius(const double& radius) { 
            m_radius = radius; 
            m_area = 2 * M_PI * (radius * m_height +  pow(radius, 2)); 
            m_volume = 2 * M_PI * radius * m_height;  
        }

        void set_height(const double& height) { 
            m_height = height; 
            m_area = 2 * M_PI * (m_radius * height +  pow(m_radius, 2)); 
            m_volume = 2 * M_PI * m_radius * height;  
        }

        double get_radius() const { return m_radius; }

        double get_height() const { return m_height; }

        void print_radius() const { std::cout << "- radius = " << get_radius() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }

        void print_height() const { std::cout << "- height = " << get_height() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }

};


class Sphere : public Shape3D {

	public:
        
        // =============================================
        // class member
        // =============================================

        double m_radius; 


        // =============================================
        // constructor and destructor
        // =============================================

		Sphere(const double& radius, const std::vector<double> center, const char* shape_um_prefix = "") : 
            m_radius{radius} {
            m_center.set_coordinates(center); 
            m_area = 4 * M_PI * pow(radius, 2); 
            m_volume = 4. * radius * pow(M_PI, 3) / 3.; 
            m_shape_um_prefix = shape_um_prefix;  
        }

        ~Sphere() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_radius(const double& radius) { 
            m_radius = radius; 
            m_area = 4 * M_PI * pow(radius, 2); 
            m_volume = 4. * radius * pow(M_PI, 3) / 3.;        
        }

        double get_radius() const { return m_radius; }

        void print_radius() const { std::cout << "- radius = " << get_radius() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }

};


class Cube : public Shape3D {

	public:
        
        // =============================================
        // class member
        // =============================================

        double m_side; 


        // =============================================
        // constructor and destructor
        // =============================================

		Cube(const double& side, const std::vector<double> center, const char* shape_um_prefix = "") : 
            m_side{side} {
            m_center = center;
            m_area = 6 * pow(side, 2);
            m_volume = pow(side, 3);
            m_shape_um_prefix = shape_um_prefix;  
        }

        ~Cube() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_side(const double& side) { 
            m_side = side; 
            m_area = 6 * pow(side, 2); 
            m_volume = pow(side, 3);        
        }

        double get_side() const { return m_side; }

        void print_side() const { std::cout << "- side = " << get_side() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }

};


class Parallelepiped : public Shape3D {

	public:
        
        // =============================================
        // class member
        // =============================================

        double m_side; 

        double m_height;

        double m_depth; 


        // =============================================
        // constructor and destructor
        // =============================================

		Parallelepiped(const double& side, const double& height, const double& depth, const std::vector<double> center, const char* shape_um_prefix = "") :
            m_side{side}, m_height{height}, m_depth{depth} {
            m_center = center;
            m_area = 2 * (side * height + side * depth + height * depth);
            m_volume = side * height * depth;
        }

        ~Parallelepiped() {}


        // =============================================
        // set and get methods
        // =============================================

        void set_side(const double& side) { 
            m_side = side; 
            m_area = 2 * (side * m_height + side * m_depth + m_height * m_depth);
            m_volume = side * m_height * m_depth;  
        }

        void set_height(const double& height) { 
            m_height = height; 
            m_area = 2 * (m_side * height + m_side * m_depth + height * m_depth);
            m_volume = m_side * height * m_depth;  
        }

        void set_depth(const double& depth) { 
            m_depth = depth; 
            m_area = 2 * (m_side * m_height + m_side * depth + m_height * depth);
            m_volume = m_side * m_height * depth;  
        }

        double get_side() const { return m_side; }

        double get_height() const { return m_height; }

        double get_depth() const { return m_depth; }

        void print_side() const { std::cout << "- side = " << get_side() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }

        void print_height() const { std::cout << "- height = " << get_height() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }
        
        void print_depth() const { std::cout << "- depth = " << get_depth() << " " << get_shape_um_prefix() << get_shape_um() << std::endl; }

};
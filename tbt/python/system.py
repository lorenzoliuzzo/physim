
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     
# last updated:    26/06/2022

class System :

    def __init__ 

template <typename T> 
class PhysicsSystemBase : Posizione {
    private:
        std::vector<T> m_objects; 
        unsigned int m_objects_count; 
        double m_kineticEnergy;
        double m_potentialEnergy;
        double m_totalEnergy; 

    public:
        PhysicsSystemBase() :
            m_kineticEnergy{},
            m_potentialEnergy{}, 
            m_totalEnergy{} {} 

        ~PhysicsSystemBase() {}

        void addObject(T x) { m_objects.push_back(x); } 

        virtual void calculateForces() const = 0;

        virtual void calculateEnergy() const = 0; 

        void setKineticEnergy(double K) { m_kineticEnergy = K; }

        void setPotentialEnergy(double U) { m_potentialEnergy = U; }

        void setTotalEnergy(double E) { m_totalEnergy = E; }

        void resetObjects() {
            m_objects.clear(); 
        }

        void resetEnergy() {
            m_kineticEnergy = 0; 
            m_potentialEnergy = 0; 
            m_totalEnergy = 0;            
        }

        int getNoObjects() const {
            return m_objects.size();
        }

        double getTotalEnergy() const {
            return m_totalEnergy;
        }

        double getPotentialEnergy() const {
            return m_potentialEnergy;
        }

        double getKineticEnergy() const {
            return m_kineticEnergy;
        }

        std::vector<T> &getObjects() {
            return m_objects;
        }

        T &getObject(unsigned pos) {
            return m_objects[pos]; 
        }

};

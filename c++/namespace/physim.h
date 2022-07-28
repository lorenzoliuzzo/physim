
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools for computational physics. 
// last updated:    27/07/2022

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip> 
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>


namespace physim {

    namespace physics {

        // namespace defining some usefull constants
        namespace constants {

            constexpr int32_t max_neg(uint32_t n_bits) { return -(int32_t(1U << (n_bits - 1))); }
            constexpr double infinity = std::numeric_limits<double>::infinity();
            constexpr double invalid_conversion = std::numeric_limits<double>::signaling_NaN();

            constexpr double pi = 3.14159265358979323846;
            constexpr double e = 2.7182818284590452353603;

        } // namespace constants

    } // namespace physics 


    namespace math {
        
        // namespace defining some usefull operation
        namespace op {

            // generate the square power of a value
            template<typename T> constexpr T square(const T& a) { return a * a; }

            // generate the cubic power of a value
            template<typename T> constexpr T cube(const T& a) { return a * a * a; }

            // generate small integer powers of a value(1, 0, -1)
            template<typename T> constexpr T pow_small(const T& val, const int& power) { 
                return (power == 1) ? val : ((power == -1) ? T(1.0) / val : T(1.0));
            }

            // generate an integer power of a number
            template<typename T> constexpr T pow(T val, const int& power) {
                return (power > 1) ? square(pow(val, power / 2)) * (power % 2 == 0 ? T(1.0) : val) :
                    (power < -1) ? T(1.0) / (square(pow(val, (-power) / 2)) * ((-power) % 2 == 0 ? T(1.0) : val)) :
                    pow_small(val, power);
            }

            // generate root power of a value
            template<typename X> X root(X value, const int& power) {
                switch (power) {
                    case 0: return X{1.0};
                    case 1: return value;
                    case -1: return X{1.0} / value;
                    case 2: 
                        if (value < X{0.0}) { return physics::constants::invalid_conversion; }
                        else return std::sqrt(value);
                    case -2:
                        if (value < X{0.0}) { return physics::constants::invalid_conversion; }
                        else return std::sqrt(X{1.0} / value);
                    case 3: return std::cbrt(value);
                    case -3: return std::cbrt(X{1.0} / value);
                    case 4: 
                        if (value < X{0.0}) { return physics::constants::invalid_conversion; }
                        else return std::sqrt(std::sqrt(value));
                    case -4: 
                        if (value < X{0.0}) { return physics::constants::invalid_conversion; }
                        else return std::sqrt(std::sqrt(X{1.0} / value));
                    default:
                        if (value < X{0.0} && power % 2 == 0) { return physics::constants::invalid_conversion; }
                        else return std::pow(value, X{1.0} / static_cast<X>(power));
                }
            }
            
            // generate the square root power of a value
            template <typename X> X sqrt(X x) { return root(x, 2); }

            // generate the cubic root power of a value
            template <typename X> X cbrt(X x) { return root(x, 3); }

            // round a value to the expected level of precision of a double
            inline double cround(double val) {
                std::uint64_t bits;
                std::memcpy(&bits, &val, sizeof(bits));
                bits += 0x800ULL;
                bits &= 0xFFFFFFFFFFFFF000ULL;
                std::memcpy(&val, &bits, sizeof(bits));
                return val;
            }

            // rounding compare for equality on double
            inline bool compare_round_equals(double val1, double val2) {
                static constexpr double half_precise_precision{5e-13};
                auto v1 = val1 - val2;
                if (v1 == 0.0 || std::fpclassify(v1) == FP_SUBNORMAL) { return true; }
                auto c1 = cround(val1);
                auto c2 = cround(val2);
                return (c1 == c2) ||
                    (cround(val2 * (1.0 + half_precise_precision)) == c1) ||
                    (cround(val2 * (1.0 - half_precise_precision)) == c1) ||
                    (cround(val1 * (1.0 + half_precise_precision)) == c2) ||
                    (cround(val1 * (1.0 - half_precise_precision)) == c2);
            }

        } // namespace op

    } // namespace math
    

    namespace physics {

        // namespace defining the units of measurements
        namespace units {
            
            // number of bits used for encoding base unit exponents 
            namespace bitwidth {

                constexpr uint32_t base_size = sizeof(uint32_t) == 8 ? 8 : 4;
                constexpr uint32_t meter{(base_size == 8) ? 8 : 4};
                constexpr uint32_t second{(base_size == 8) ? 8 : 4};
                constexpr uint32_t kilogram{(base_size == 8) ? 6 : 3};
                constexpr uint32_t ampere{(base_size == 8) ? 6 : 3};
                constexpr uint32_t candela{(base_size == 8) ? 4 : 2};
                constexpr uint32_t kelvin{(base_size == 8) ? 6 : 3};
                constexpr uint32_t mole{(base_size == 8) ? 4 : 2};

            } // namespace bitwith


            // class representing the seven SI base units 
            class unit_data {

                private: 

                    signed int meter_ : bitwidth::meter;
                    signed int second_ : bitwidth::second;  
                    signed int kilogram_ : bitwidth::kilogram;
                    signed int ampere_ : bitwidth::ampere;
                    signed int candela_ : bitwidth::candela;  
                    signed int kelvin_ : bitwidth::kelvin;
                    signed int mole_ : bitwidth::mole;


                public:

                    // the seven SI base units 
                    enum base {
                        Meter = 0,
                        Second = 1,
                        Kilogram = 2,
                        Ampere = 3,
                        Kelvin = 4,
                        Mole = 5,
                        Candela = 6
                    };

                    static constexpr uint32_t bits[7] = { 
                        bitwidth::meter, bitwidth::second, 
                        bitwidth::kilogram, bitwidth::ampere, 
                        bitwidth::kelvin, bitwidth::mole, 
                        bitwidth::candela
                    };  


                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // constructor from powers
                    constexpr unit_data(const int& meters, const int& seconds, const int& kilograms, 
                                        const int& amperes, const int& kelvins,  const int& moles, const int& candelas) :
                        meter_(meters), second_(seconds), kilogram_(kilograms), ampere_(amperes),
                        kelvin_(kelvins), mole_(moles), candela_(candelas) {};

                    explicit constexpr unit_data(std::nullptr_t) :
                        meter_(constants::max_neg(bitwidth::meter)), second_(constants::max_neg(bitwidth::second)), 
                        kilogram_(constants::max_neg(bitwidth::kilogram)), ampere_(constants::max_neg(bitwidth::ampere)),
                        kelvin_(constants::max_neg(bitwidth::kelvin)), mole_(constants::max_neg(bitwidth::mole)), 
                        candela_(constants::max_neg(bitwidth::candela)) {}


                    // =============================================
                    // operators
                    // ============================================= 

                    // perform a multiply operation by adding the powers together
                    constexpr unit_data operator*(const unit_data& other) const {
                        return { 
                            meter_ + other.meter_,
                            second_ + other.second_,
                            kilogram_ + other.kilogram_,
                            ampere_ + other.ampere_,
                            kelvin_ + other.kelvin_,
                            mole_ + other.mole_,
                            candela_ + other.candela_,
                        };
                    }

                    // perform a division operation by subtract the powers together
                    constexpr unit_data operator/(const unit_data& other) const {
                        return { 
                            meter_ - other.meter_,
                            second_ - other.second_,
                            kilogram_ - other.kilogram_,
                            ampere_ - other.ampere_,
                            kelvin_ - other.kelvin_,
                            mole_ - other.mole_,
                            candela_ - other.candela_,
                        };
                    }

                    // invert the unit
                    constexpr unit_data inv() const {
                        return { 
                            -meter_,
                            -second_,
                            -kilogram_,
                            -ampere_,
                            -kelvin_,
                            -mole_,
                            -candela_
                        };
                    }

                    // take a unit_data to some power
                    constexpr unit_data pow(const int& power) const { 
                        return { 
                            meter_ * power,
                            (second_ * power) + root_Hertz_modifier(power),
                            kilogram_ * power,
                            ampere_ * power,
                            kelvin_ * power,
                            mole_ * power,
                            candela_ * power
                        };
                    }
                    
                    // take some root of a unit_data
                    constexpr unit_data root(const int& power) const {
                        if (has_valid_root(power)) return unit_data( meter_ / power,
                                                                    second_ / power,
                                                                    kilogram_ / power,
                                                                    ampere_ / power,
                                                                    kelvin_ / power,
                                                                    mole_ / power,
                                                                    candela_ / power);
                        else exit(-11);
                    }
                    
                    
                    // =============================================
                    // check methods
                    // ============================================= 

                    // comparison operators
                    constexpr bool operator==(const unit_data& other) const {
                        return meter_ == other.meter_ && second_ == other.second_ &&
                            kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                            candela_ == other.candela_ && kelvin_ == other.kelvin_ && mole_ == other.mole_;         
                    }

                    constexpr bool operator!=(const unit_data& other) const {
                        return !(*this == other);
                    }

                    // check if the units have the same base unit 
                    constexpr bool has_same_base(const unit_data& other) const {
                        return *this == other;
                    }

                    // check if the unit is empty
                    constexpr bool empty() const {
                        return meter_ == 0 && second_ == 0 && kilogram_ == 0 &&
                            ampere_ == 0 && candela_ == 0 && kelvin_ == 0 && mole_ == 0;
                    }
                    

                    // =============================================
                    // get methods
                    // =============================================
                    
                    constexpr int meter() const { return meter_; }

                    constexpr int second() const { return second_; }

                    constexpr int kg() const { return kilogram_; }
                                        
                    constexpr int ampere() const { return ampere_; }
                    
                    constexpr int kelvin() const { return kelvin_; }
                    
                    constexpr int mole() const { return mole_; }
                    
                    constexpr int candela() const { return candela_; }

                    constexpr int unit_type_count() const {
                        return ((meter_ != 0) ? 1 : 0) + ((second_ != 0) ? 1 : 0) +
                            ((kilogram_ != 0) ? 1 : 0) + ((ampere_ != 0) ? 1 : 0) +
                            ((kelvin_ != 0) ? 1 : 0) + ((mole_ != 0) ? 1 : 0) + 
                            ((candela_ != 0) ? 1 : 0);
                    }

                    constexpr unit_data base_units() const { return *this; }

                    constexpr void print() const {
                        if (meter_ != 0 && meter_ != 1) std::cout << "m^" << meter_ << "\n"; 
                        if (meter_ == 1) std::cout << "m" << "\n";
                        if (second_ != 0 && second_ != 1) std::cout << "s^" << second_ << "\n"; 
                        if (second_ == 1) std::cout << "s" << "\n"; 
                        if (kilogram_ != 0 && kilogram_ != 1) std::cout << "kg^" << kilogram_ << "\n"; 
                        if (kilogram_ == 1) std::cout << "kg" << "\n"; 
                        if (ampere_ != 0 && ampere_ != 1) std::cout << "A^" << ampere_ << "\n"; 
                        if (ampere_ == 1) std::cout << "A" << "\n"; 
                        if (kelvin_ != 0 && kelvin_ != 1) std::cout << "K^" << kelvin_ << "\n"; 
                        if (kelvin_ == 1) std::cout << "K" << "\n";
                        if (mole_ != 0 && mole_ != 1) std::cout << "mol^" << mole_ << "\n"; 
                        if (mole_ == 1) std::cout << "mol"  << "\n"; 
                        if (candela_ != 0 && candela_ != 1) std::cout << "cd^" << candela_ << "\n"; 
                        if (candela_ == 1) std::cout << "cd" << "\n"; 
                    }


                private: 

                    // check if the base_unit has a valid root
                    constexpr bool has_valid_root(const int& power) const {
                        return meter_ % power == 0 && second_ % power == 0 &&
                            kilogram_ % power == 0 && ampere_ % power == 0 &&
                            candela_ % power == 0 && kelvin_ % power == 0 &&
                            mole_ % power == 0;
                    }      
                    
                    // to handle a few weird operations that operate on square_root Hz
                    constexpr int root_Hertz_modifier(const int& power) const {
                        return (second_ * power == 0 || 
                                power % 2 != 0) ? 0 : (power / 2) * ((second_ < 0) || 
                                (power < 0) ? 9 : -9);
                    }    

            }; // class unit_base


            // class defining a basic unit module 
            class unit_prefix {

                protected: 

                    // =============================================
                    // class members
                    // ============================================= 
                    
                    double multiplier_{1.0};  

                    const char prefix_ = '-'; 


                public: 

                    // =============================================
                    // constructors
                    // ============================================= 

                    // default constructor
                    constexpr unit_prefix() noexcept {};

                    explicit constexpr unit_prefix(const double& mult) noexcept : 
                        multiplier_{mult} {} 

                    explicit constexpr unit_prefix(const double& mult, const char& pref) noexcept : 
                        multiplier_{mult}, prefix_(pref) {} 

                    explicit constexpr unit_prefix(const char& pref, const double& mult) noexcept : 
                        multiplier_{mult}, prefix_(pref) {} 
                        

                    // =============================================
                    // operators
                    // ============================================= 

                    // perform a multiply operation by adding the powers together
                    constexpr unit_prefix operator*(const unit_prefix& other) const {
                        return unit_prefix(multiplier_ + other.multiplier(), prefix_ + other.prefix_);
                    }

                    // perform a division operation by subtract the powers together
                    constexpr unit_prefix operator/(const unit_prefix& other) const {
                        return unit_prefix(multiplier_ / other.multiplier(), prefix_);
                    }

                    // invert the unit
                    constexpr unit_prefix inv() const {
                        return unit_prefix(1 / multiplier_, prefix_);
                    }

                    // take a unit_prefix to some power
                    constexpr unit_prefix pow(const int& power) const { 
                        return unit_prefix(math::op::pow(multiplier_, power), prefix_);
                    }
                    
                    // take some root of a unit_prefix
                    unit_prefix root(const int& power) const {
                        return unit_prefix(math::op::root(multiplier_, power), prefix_); 
                    }
                    
                    
                    // =============================================
                    // check methods
                    // ============================================= 

                    // comparison operators
                    constexpr bool operator==(const unit_prefix& other) const {
                        if (multiplier_ == other.multiplier_) return true;    
                        else return math::op::compare_round_equals(multiplier_, other.multiplier_);
                    }

                    constexpr bool operator!=(const unit_prefix& other) const {
                        return !(*this == other);
                    }

                    // check if the units have the same prefix unit 
                    constexpr bool has_same_prefix(const unit_prefix& other) const {
                        return *this == other;
                    }

                    // check if the multiplier is nan
                    inline bool is_nan(const unit_prefix& pref) {
                        return std::isnan(pref.multiplier_);
                    }

                    // checks that the multiplier is finite
                    inline bool is_finite(const unit_prefix& pref) {
                        return std::isfinite(pref.multiplier_);
                    }

                    // check if the multiplier is infinite
                    inline bool is_inf(const unit_prefix& pref) {
                        return std::isinf(pref.multiplier_);
                    }


                    // =============================================
                    // get methods
                    // ============================================= 
                    
                    constexpr double multiplier() const { return multiplier_; }

                    constexpr unit_prefix prefix() const { return *this; }

                    constexpr void print() const { 
                        if (prefix_ != '-') std::cout << prefix_; 
                        else {
                            if (multiplier_ == 1.0) std::cout << ""; 
                            if (multiplier_ == 1e-1) std::cout << "d"; 
                            if (multiplier_ == 1e-2) std::cout << "c"; 
                            if (multiplier_ == 1e-3) std::cout << "m"; 
                            if (multiplier_ == 1e-6) std::cout << "u"; 
                            if (multiplier_ == 1e-9) std::cout << "n"; 
                            if (multiplier_ == 1e-12) std::cout << "p"; 
                            if (multiplier_ == 1e-15) std::cout << "f"; 
                            if (multiplier_ == 1e-18) std::cout << "a"; 
                            if (multiplier_ == 1e-21) std::cout << "z"; 
                            if (multiplier_ == 1e-24) std::cout << "y"; 
                            if (multiplier_ == 1e2) std::cout << "h"; 
                            if (multiplier_ == 1e3) std::cout << "k"; 
                            if (multiplier_ == 1e6) std::cout << "M"; 
                            if (multiplier_ == 1e9) std::cout << "G"; 
                            if (multiplier_ == 1e12) std::cout << "T"; 
                            if (multiplier_ == 1e15) std::cout << "P"; 
                            if (multiplier_ == 1e18) std::cout << "E"; 
                            if (multiplier_ == 1e21) std::cout << "Z"; 
                            if (multiplier_ == 1e24) std::cout << "Y"; 
                        }
                    }

            }; // class unit_prefix


            // class defining a basic unit module 
            class unit : public unit_data, public unit_prefix {

                private:

                    // =============================================
                    // class member 
                    // ============================================= 

                    std::string unit_name = "-";  


                public:

                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // default constructor
                    unit() noexcept : 
                        unit_data(0, 0, 0, 0, 0, 0, 0), unit_prefix() {}

                    // constructor from unit_data 
                    explicit unit(const unit_data& unit_base) noexcept : 
                        unit_data(unit_base), unit_prefix() {}

                    // constructor from unit_data and unit_prefix
                    explicit unit(const unit_data& unit_base, const unit_prefix& prefix) noexcept : 
                        unit_data(unit_base), unit_prefix(prefix) {}

                    // constructor from unit_prefix and unit_data
                    explicit unit(const unit_prefix& prefix, const unit_data& unit_base) noexcept : 
                        unit_data(unit_base), unit_prefix(prefix) {}

                    // constructor from unit_data and a double
                    unit(const unit_data& unit_base, const double& mult) noexcept : 
                        unit_data(unit_base), unit_prefix(mult) {}

                    // constructor from double with a unit_data
                    unit(const double& mult, const unit_data& unit_base) noexcept :
                        unit_data(unit_base), unit_prefix(mult) {}

                    // constructor from unit_data and std::string name
                    unit(const unit_data& unit_base, std::string name) noexcept : 
                        unit_data(unit_base), unit_prefix(), unit_name{name} {}

                    // constructor from std::string name and unit_data
                    unit(std::string name, const unit_data& unit_base) noexcept : 
                        unit_data(unit_base), unit_prefix(), unit_name{name} {}
                        
                    // constructor from unit_data, a unit_prefix and a std::string name
                    unit(const unit_data& unit_base, const unit_prefix& prefix, std::string name) noexcept : 
                        unit_data(unit_base), unit_prefix(prefix), unit_name{name} {}

                    // constructor from a unit_prefix, an unit_data and a std::string name
                    unit(const unit_prefix& prefix, const unit_data& unit_base, std::string name) noexcept :
                        unit_data(unit_base), unit_prefix(prefix), unit_name{name} {}

                    // constructor from unit_data, a double and a std::string name
                    unit(const unit_data& unit_base, const double& mult, std::string name) noexcept : 
                        unit_data(unit_base), unit_prefix(mult), unit_name{name} {}

                    // constructor from a double, an unit_data and a std::string name
                    unit(const double& mult, const unit_data& unit_base, std::string name) noexcept :
                        unit_data(unit_base), unit_prefix(mult), unit_name{name} {}

                    // constructor from unit
                    unit(const unit& other) noexcept :
                        unit(other.base_units(), other.prefix(), other.unit_string()) {} 

                    // constructor from unit and unit_prefix
                    unit(const unit& other, const unit_prefix& prefix) noexcept : 
                        unit(other.base_units(), prefix, other.unit_string()) {} 

                    // constructor from unit_prefix and unit
                    unit(const unit_prefix& prefix, const unit& other) noexcept : 
                        unit(other.base_units(), prefix, other.unit_string()) {} 

                    // constructor from unit and a double
                    unit(const unit& other, const double& mult) noexcept : 
                        unit(other.base_units(), unit_prefix(mult), other.unit_string()) {}

                    // constructor from double with a unit
                    unit(const double& mult, const unit& other) noexcept :
                        unit(other.base_units(), unit_prefix(mult), other.unit_string()) {}

                    // constructor from unit and std::string name
                    unit(const unit& other, std::string name) noexcept : 
                        unit(other.base_units(), other.prefix(), name) {}

                    // constructor from std::string name and unit
                    unit(std::string name, const unit& other) noexcept : 
                        unit(other.base_units(), other.prefix(), name) {}
                        
                    // constructor from unit, a unit_prefix and a std::string name
                    unit(const unit& other, const unit_prefix& prefix, std::string name) noexcept : 
                        unit(other.base_units(), unit_prefix(prefix), name) {}

                    // constructor from a double, an unit and a std::string name
                    unit(const unit_prefix& prefix, const unit& other, std::string name) noexcept :
                        unit(other.base_units(), unit_prefix(prefix), name) {}                    
                    
                    // constructor from unit, a double and a std::string name
                    unit(const unit& other, const double& mult, std::string name) noexcept : 
                        unit(other.base_units(), unit_prefix(mult), name) {}

                    // constructor from a double, an unit and a std::string name
                    unit(const double& mult, const unit& other, std::string name) noexcept :
                        unit(other.base_units(), unit_prefix(mult), name) {}


                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a unit
                    unit inv() const { return unit(unit_data::inv(), unit_prefix::inv()); }

                    // multiply with another unit
                    unit operator*(const unit& other) const {
                        if (unit_name != "-") return unit(unit_data::operator*(other.base_units()), unit_prefix::operator*(other.prefix()), unit_name + other.unit_string());
                        else return unit(unit_data::operator*(other.base_units()), unit_prefix::operator*(other.prefix()));
                    }

                    // division operator
                    unit operator/(const unit& other) const {
                        if (unit_name != "-") return unit(unit_data::operator*(other.base_units()), unit_prefix::operator/(other.prefix()), unit_name + "/" + other.unit_string());
                        else return unit(unit_data::operator/(other.base_units()), unit_prefix::operator/(other.prefix()));
                    }

                    // take a unit to a power
                    unit pow(const int& power) const {
                        if (unit_name != "-") return unit(unit_data::pow(power), unit_prefix::pow(power), unit_name + "^" + std::to_string(power));
                        else return unit(unit_data::pow(power), unit_prefix::pow(power));
                    }

                    // take a unit to a root power
                    unit root(const int& power) const {
                        if (unit_name != "-") {
                            return unit(unit_data::root(power), unit_prefix::root(power), unit_name + "^(1/" + std::to_string(power) + ")");
                        }
                        else return unit(unit_data::root(power), unit_prefix::root(power));                        }


                    // =============================================
                    // checks
                    // ============================================= 

                    // equality operator
                    bool operator==(const unit& other) const {
                        if (unit_data::operator!=(other.base_units())) { return false; }
                        if (unit_prefix::operator==(other.prefix())) { return true; }
                        return math::op::compare_round_equals(multiplier(), other.multiplier());
                    }

                    // equality operator
                    bool operator!=(const unit& other) const { return !operator==(other); }

                    // test for exact numerical equivalence
                    bool is_exactly_the_same(const unit& other) const {
                        return unit_data::operator==(other.base_units()) && unit_prefix::operator==(other.prefix());
                    }

                    
                    // =============================================
                    // get methods
                    // ============================================= 

                    // get the unit_name
                    inline std::string unit_string() const { return unit_name; }

                    // get the unit
                    inline unit as_unit() const { return *this; }
                    
                    // get a rounded value of the multiplier rounded to the defined precision
                    double cround() const { return math::op::cround(multiplier()); }

                    void print() const {
                        if (unit_name != "-") { std::cout << unit_name << "\n"; }
                        else {
                            unit_prefix::print(); 
                            unit_data::print();
                        }  
                    }

            }; // class unit     


            // units declarations
            namespace defined {
                
                // SI_base
                namespace SI_base {

                    constexpr unit_data meter(1, 0, 0, 0, 0, 0, 0);
                    constexpr unit_data second(0, 1, 0, 0, 0, 0, 0);
                    constexpr unit_data kilogram(0, 0, 1, 0, 0, 0, 0);
                    constexpr unit_data Ampere(0, 0, 0, 1, 0, 0, 0);
                    constexpr unit_data Kelvin(0, 0, 0, 0, 1, 0, 0);
                    constexpr unit_data mol(0, 0, 0, 0, 0, 1, 0);
                    constexpr unit_data candela(0, 0, 0, 0, 0, 0, 1);

                } // namespace SI

                // SI_prefix
                namespace SI_prefix {

                    constexpr unit_prefix deci(1e-1, 'd'); 
                    constexpr unit_prefix centi(1e-2, 'c'); 
                    constexpr unit_prefix milli(1e-3, 'm'); 
                    constexpr unit_prefix micro(1e-6, 'u'); 
                    constexpr unit_prefix nano(1e-9, 'n'); 
                    constexpr unit_prefix pico(1e-12, 'p'); 
                    constexpr unit_prefix femto(1e-15, 'f'); 
                    constexpr unit_prefix atto(1e-18, 'a'); 
                    constexpr unit_prefix zepto(1e-21, 'z'); 
                    constexpr unit_prefix yocto(1e-24, 'y'); 
                    constexpr unit_prefix hecto(1e2, 'h'); 
                    constexpr unit_prefix kilo(1e3, 'k'); 
                    constexpr unit_prefix mega(1e6, 'M'); 
                    constexpr unit_prefix giga(1e9, 'G'); 
                    constexpr unit_prefix tera(1e12, 'T'); 
                    constexpr unit_prefix peta(1e15, 'P'); 
                    constexpr unit_prefix exa(1e18, 'E'); 
                    constexpr unit_prefix zetta(1e21, 'Z'); 
                    constexpr unit_prefix yotta(1e24, 'Y'); 

                } // namespace SI_prefix
    
                // SI 
                unit m = unit(SI_base::meter);
                unit s = unit(SI_base::second, "s");
                unit kg = unit(SI_base::kilogram, "kg");
                unit A = unit(SI_base::Ampere, "A");
                unit K = unit(SI_base::Kelvin, "K");
                unit mol = unit(SI_base::mol, "mol"); 
                unit cd = unit(SI_base::candela, "cd");

                // some unitless numbers
                unit one;
                unit hundred = unit(100.0, one);
                unit ten = unit(10.0, one);
                unit percent(one, 0.01, "%");
                unit infinite(unit_data(0, 0, 0, 0, 0, 0, 0), constants::infinity, "inf");
                unit neginfinite(unit_data(0, 0, 0, 0, 0, 0, 0), -constants::infinity, "-inf");
                unit nan(unit_data(0, 0, 0, 0, 0, 0, 0), constants::invalid_conversion, "NaN");

                // some specialized units
                unit defunit(unit_data(0, 0, 0, 0, 0, 0, 0));
                unit invalid(unit_data(nullptr), constants::invalid_conversion, "NaN");
                unit error(unit_data(nullptr));
            
                // SI_derived units
                unit hertz(unit_data(0, 0, -1, 0, 0, 0, 0), "Hz");
                unit Hz = hertz;
                unit volt(unit_data(2, 1, -3, -1, 0, 0, 0), "V");
                unit V = volt;
                unit newton(unit_data(1, 1, -2, 0, 0, 0, 0), "N");
                unit N = newton;
                unit Pa(unit_data(-1, 1, -2, 0, 0, 0, 0), "Pa");
                unit pascal = Pa;
                unit joule(unit_data(2, 1, -2, 0, 0, 0, 0), "J");
                unit J = joule;
                unit watt(unit_data(2, 1, -3, 0, 0, 0, 0), "W");
                unit W = watt;
                unit coulomb(unit_data(0, 0, 1, 1, 0, 0, 0), "C");
                unit C = coulomb;
                unit farad(unit_data(-2, -1, 4, 2, 0, 0, 0), "F");
                unit F = farad;
                unit weber(unit_data(2, 1, -2, -1, 0, 0, 0), "Wb");
                unit Wb = weber;
                unit tesla(unit_data(0, 1, -2, -1, 0, 0, 0), "T");
                unit T = tesla;
                unit henry(unit_data(2, 1, -2, -2, 0, 0, 0), "H");                    
                unit H = henry;
                unit mps(m / s);
                unit mpss(m / s.pow(2)); 

                // distance units
                unit km(SI_prefix::kilo, m);
                unit dm(SI_prefix::deci, m);
                unit cm(SI_prefix::centi, m);
                unit mm(SI_prefix::milli, m);
                unit um(SI_prefix::micro, m);
                unit nm(SI_prefix::nano, m);

                // time units
                unit ms(SI_prefix::milli, s);
                unit us(SI_prefix::micro, s);
                unit ns(SI_prefix::nano, s);
                unit min(60.0, s);
                unit hr(60.0, min, "hr");
                unit day(24.0, hr, "day");
                unit yr(8760.0, hr, "yr");  // median calendar year;
                unit sday(365.24 / 366.24, day, "sday");  // sidereal day
                unit syr(365.256363004, day, "syr");  // sidereal year

                // mass units
                unit g(SI_prefix::milli, kg, "g");
                unit mg(SI_prefix::micro, kg, "mg");

                // volume units
                unit L{1, dm.pow(3), "L"};
                unit dL{0.1, L, "dL"};
                unit cL{0.01, L, "cL"};
                unit mL{0.001, L, "mL"};                    

            } // namespace defined


        } // namespace units

    } // namespace physics


    namespace math {
        
        // namespace defining some usefull operation
        namespace op {

            // generate a unit which is an integer power of another
            inline physics::units::unit pow(const physics::units::unit& u, const int& power) { return u.pow(power); }

            // generate the square of an physics::units::unit
            inline physics::units::unit square(const physics::units::unit& u) { return u.pow(2); }

            // generate the cube of an physics::units::unit
            inline physics::units::unit cube(const physics::units::unit& u) { return u.pow(3); }

            // generate the root of an physics::units::unit
            physics::units::unit root(const physics::units::unit& un, const int& power) {
                if (power == 0) { return physics::units::defined::one; }
                if (un.multiplier() < 0.0 && power % 2 == 0) { return physics::units::defined::invalid; }
                return physics::units::unit{ un.base_units().root(power), math::op::root(un.multiplier(), power) };
            }

            // generate the square root of an physics::units::unit
            inline physics::units::unit sqrt(const physics::units::unit& u) { return root(u, 2); }

            // generate the cubic root of an physics::units::unit
            inline physics::units::unit cbrt(const physics::units::unit& u) { return root(u, 3); }

            // generate a conversion factor between two physics::units::units in a constexpr function, the
            // physics::units::units will only convert if they have the same base physics::units::unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(UX start, UX2 result) { return quick_convert(1.0, start, result); }

            // generate a conversion factor between two physics::units::units in a constexpr function, the physics::units::units will only convert if they have the same base physics::units::unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(double val, const UX& start, const UX2& result) {
                static_assert(std::is_same<UX, physics::units::unit>::value || std::is_same<UX, physics::units::unit>::value,
                    "convert argument types must be physics::units::unit or physics::units::unit");
                static_assert(std::is_same<UX2, physics::units::unit>::value || std::is_same<UX2, physics::units::unit>::value,
                    "convert argument types must be physics::units::unit or physics::units::unit");
                return (start.base_units() == result.base_units()) ?  val * start.multiplier() / result.multiplier() : 
                    physics::constants::invalid_conversion;
            }

            // convert a value from one physics::units::unit base to another
            template<typename UX, typename UX2>
            double convert(double val, const UX& start, const UX2& result) {
                static_assert(std::is_same<UX, physics::units::unit>::value || std::is_same<UX, physics::units::unit>::value,
                    "convert argument types must be physics::units::unit or physics::units::unit");
                static_assert(std::is_same<UX2, physics::units::unit>::value || std::is_same<UX2, physics::units::unit>::value, 
                    "convert argument types must be physics::units::unit or physics::units::unit");     
                if (start == result) { return val; }
                if (start.base_units() == result.base_units()) { return val * start.multiplier() / result.multiplier(); }
                auto base_start = start.base_units();
                auto base_result = result.base_units();
                if (base_start.has_same_base(base_result)) { return val * start.multiplier() / result.multiplier(); }
                if (base_start.has_same_base(base_result.inv())) { return 1.0 / (val * start.multiplier() * result.multiplier()); }
                return physics::constants::invalid_conversion;
            }     
                            
            // generate a conversion factor between two physics::units::units
            template<typename UX, typename UX2>
            double convert(const UX& start, const UX2& result) { return convert(1.0, start, result); }

        } // namespace op

    } // namespace math
    

    namespace physics {

        // namespace defining some usefull measurement structures
        namespace measurements {

            // class using units and double precision
            class measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0};

                    units::unit units_;


                public:

                    // =============================================                                                                                         
                    // constructors
                    // =============================================  

                    // default constructor
                    measurement() noexcept {};

                    // constructor from a value and a unit 
                    measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    // constructor from a measurement
                    measurement(const measurement& other) :
                        value_(other.value()), units_(other.units()) {}

                    // constructor from a measurement
                    measurement(measurement&& val) noexcept :
                        value_(val.value()), units_(val.units()) {}
                        
                    // destructor
                    ~measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  
        
                    measurement operator*(const measurement& other) const {
                        return { value_ * other.value_, units_ * other.units_ };
                    }
                    
                    measurement operator*(double val) const {
                        return { value_ * val, units_ };
                    }
                    
                    measurement operator/(const measurement& other) const {
                        return { value_ / other.value_, units_ / other.units_ };
                    }

                    measurement operator/(double val) const {
                        return { value_ / val, units_ };
                    }
                    
                    constexpr bool operator==(double val) const {
                        return (value_ == val) ? true : math::op::compare_round_equals(value_, val);
                    }

                    constexpr bool operator!=(double val) const { return !operator==(val); }

                    constexpr bool operator>(double val) const { return value_ > val; }

                    constexpr bool operator<(double val) const { return value_ < val; }

                    constexpr bool operator>=(double val) const {
                        return (value_ >= val) ? true : operator==(val);
                    }

                    constexpr bool operator<=(double val) const {
                        return value_ <= val ? true : operator==(val);
                    }

                    measurement& operator=(const measurement& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }
                
                    measurement& operator=(measurement&& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    measurement& operator=(double val) noexcept {
                        value_ = val;
                        return *this;
                    }            

                    bool operator==(const measurement& other) const {
                        return value_equality_check((units_ == other.units()) ? other.value() : other.value_as(units_));
                    }

                    bool operator!=(const measurement& other) const {
                        return !value_equality_check((units_ == other.units()) ? other.value() : other.value_as(units_));
                    }

                    bool operator>(const measurement& other) const {
                        return value_ > other.value_as(units_);
                    }

                    bool operator<(const measurement& other) const {
                        return value_ < other.value_as(units_);
                    }

                    bool operator>=(const measurement& other) const {
                        double val = other.value_as(units_);
                        return (value_ > val) ? true : value_equality_check(val);
                    }

                    bool operator<=(const measurement& other) const {
                        double val = other.value_as(units_);
                        return (value_ < val) ? true : value_equality_check(val);
                    }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return { math::op::convert(value_, units_, newUnits), newUnits };
                    }

                    // convert a unit into its base units
                    measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::unit(units_.base_units()) };
                    }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? value_ : math::op::convert(value_, units_, desired_units);
                    }

                    // get the numerical component of the measurement
                    constexpr double value() const { return value_; }                        
                    
                    // set the numerical component of the measurement
                    constexpr void value(const double& value)  { value_ = value; }

                    // get the unit component of a measurement
                    units::unit units() const { return units_; }

                    // convert the measurement to a single unit
                    units::unit as_unit() const { return {value_, units_ }; }

                    // print the measurement
                    inline void print() const { 
                        std::cout << value() << " "; 
                        units_.print(); 
                    }

                private:

                    // does a numerical equality check on the value accounting for tolerances
                    bool value_equality_check(double otherval) const {
                        return (value_ == otherval) ?
                            true : math::op::compare_round_equals(value_, otherval);
                    } 

            }; // class measurement


            // class using a fixed unit and a value
            class fixed_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0};  

                    const units::unit units_;  


                public:

                    // =============================================                                                                                         
                    // constructors
                    // =============================================  

                    // constructor from a value and a unit 
                    fixed_measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    // constructor from a measurement
                    explicit fixed_measurement(const measurement& val) :
                        value_(val.value()), units_(val.units()) {}

                    // constructor from a fixed measurement
                    fixed_measurement(const fixed_measurement& val) :
                        value_(val.value()), units_(val.units()) {}

                    // constructor from a fixed measurement
                    fixed_measurement(fixed_measurement&& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    // destructor
                    ~fixed_measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    measurement operator*(const measurement& other) const {
                        return {value_ * other.value(), units_ * other.units()};
                    }
                    
                    fixed_measurement operator*(double val) const {
                        return {value_ * val, units_ };
                    }

                    measurement operator/(const measurement& other) const {
                        return {value_ / other.value(), units_ / other.units()};
                    }

                    fixed_measurement operator/(double val) const {
                        return {value_ / val, units_ };
                    }

                    fixed_measurement operator+(const measurement& other) const {
                        return {value_ + other.value_as(units_), units_ };
                    }

                    fixed_measurement operator-(const measurement& other) const {
                        return {value_ - other.value_as(units_), units_ };
                    }

                    fixed_measurement operator+(double val) const {
                        return {value_ + val, units_ };
                    }

                    fixed_measurement operator-(double val) const {
                        return {value_ - val, units_ };
                    }

                    fixed_measurement& operator+=(double val) {
                        value_ += val;
                        return *this;
                    }

                    fixed_measurement& operator-=(double val) {
                        value_ -= val;
                        return *this;
                    }

                    fixed_measurement& operator*=(double val) {
                        value_ *= val;
                        return *this;
                    }

                    fixed_measurement& operator/=(double val) {
                        value_ /= val;
                        return *this;
                    }

                    fixed_measurement& operator=(const measurement& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    fixed_measurement& operator=(const fixed_measurement& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    fixed_measurement& operator=(fixed_measurement&& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    fixed_measurement& operator=(double val) noexcept {
                        value_ = val;
                        return *this;
                    }

                    constexpr bool operator==(double val) const {
                        return (value_ == val) ? true : math::op::compare_round_equals(value_, val);
                    }

                    constexpr bool operator!=(double val) const { return !operator==(val); }

                    constexpr bool operator>(double val) const { return value_ > val; }

                    constexpr bool operator<(double val) const { return value_ < val; }

                    constexpr bool operator>=(double val) const {
                        return (value_ >= val) ? true : operator==(val);
                    }

                    constexpr bool operator<=(double val) const {
                        return value_ <= val ? true : operator==(val);
                    }

                    bool operator==(const fixed_measurement& val) const {
                        return operator==((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator!=(const fixed_measurement& val) const {
                        return operator!=((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator==(const measurement& val) const {
                        return operator==((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator!=(const measurement& val) const {
                        return operator!=((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator>(const measurement& val) const {
                        return operator>((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator<(const measurement& val) const {
                        return operator<((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator>=(const measurement& val) const {
                        return operator>=((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator<=(const measurement& val) const {
                        return operator<=((units_ == val.units()) ? val.value() : val.value_as(units_));
                    }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // direct conversion operator
                    operator measurement() { return {value_, units_ }; }

                    // set the numerical value
                    constexpr void value(double val) { value_ = val; }

                    // get the numerical value
                    constexpr double value() const { return value_; }

                    // get the unit
                    units::unit units() const { return units_; }

                    // convert the measurement to a unit
                    units::unit as_unit() const { return {value_, units_ }; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? value_ : math::op::convert(value_, units_, desired_units);
                    }

                    // convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return { math::op::convert(value_, units_, newUnits), newUnits };
                    }

                    // convert a unit into its base units
                    measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::unit(units_.base_units()) };
                    }

                    // print the measurement
                    inline void print() const { 
                        std::cout << value() << " "; 
                        units_.print(); 
                    }

            }; // class fixed_measurement
    

            // class using fixed units, a value and an uncertain value
            class uncertain_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0}; 

                    double uncertainty_{0.0};

                    units::unit units_;       
                

                public: 

                    // =============================================                                                                                         
                    // constructors
                    // =============================================  

                    // default constructor
                    uncertain_measurement() = default;
                    
                    // constructor from a value, uncertainty, and unit
                    uncertain_measurement(double val, double uncertainty_val, const units::unit& base) noexcept :
                        value_(val), uncertainty_(uncertainty_val), units_(base) {}

                    // constructpr from a value and an unit, assuming the uncertainty is 0
                    explicit uncertain_measurement(double val, const units::unit& base) noexcept :
                        value_(val), units_(base) {}

                    // constructor from a regular measurement and uncertainty value
                    explicit uncertain_measurement(const measurement& val, double uncertainty_val) noexcept : 
                        value_(val.value()), uncertainty_(uncertainty_val), units_(val.units()) {}

                    // constructor from a regular measurement and an uncertainty measurement
                    explicit uncertain_measurement(const measurement& val, const measurement& uncertainty_meas) noexcept :
                        value_(val.value()), uncertainty_(uncertainty_meas.value_as(val.units())), units_(val.units()) {}


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    // compute a product and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator*(const uncertain_measurement& other) const {
                        double tval1 = uncertainty_ / value_;
                        double tval2 = other.uncertainty_ / other.value_;
                        double ntol = std::sqrt(tval1 * tval1 + tval2 * tval2);
                        double nval = value_ * other.value_;
                        return { nval, nval * ntol, units_ * other.units() };
                    }

                    // perform a multiplication with uncertain measurements using the simple method for uncertainty propagation
                    uncertain_measurement simple_product(const uncertain_measurement& other) const {
                        double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                        double nval = value_ * other.value_;
                        return { nval, nval * ntol, units_ * other.units() };
                    }

                    // multiply with another measurement equivalent to uncertain_measurement multiplication with 0 uncertainty
                    uncertain_measurement operator*(const measurement& other) const {
                        return { 
                            value() * other.value(),
                            other.value() * uncertainty(),
                            units_ * other.units() };
                    }

                    uncertain_measurement operator*(const units::unit& other) const {
                        return { value_, uncertainty_, units_ * other};
                    }

                    uncertain_measurement operator*(double val) const {
                        return { value_ * val, uncertainty_ * val, units_ };
                    }

                    // compute a unit division and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator/(const uncertain_measurement& other) const {
                        double tval1 = uncertainty_ / value_;
                        double tval2 = other.uncertainty_ / other.value_;
                        double ntol = std::sqrt(tval1 * tval1 + tval2 * tval2);
                        double nval = value_ / other.value_;
                        return { nval, nval * ntol, units_ / other.units() };
                    }

                    // division operator propagate uncertainty using simple method
                    uncertain_measurement simple_divide(const uncertain_measurement& other) const {
                        double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                        double nval = value_ / other.value_;
                        return { nval, nval * ntol, units_ / other.units() };
                    }

                    uncertain_measurement operator/(const measurement& other) const {
                        return { static_cast<double>(value() / other.value()),
                            static_cast<double>(uncertainty() / other.value()), units_ / other.units() };
                    }

                    uncertain_measurement operator/(const units::unit& other) const {
                        return { value_, uncertainty_, units_ / other};
                    }

                    uncertain_measurement operator/(double val) const {
                        return { value_ / val, uncertainty_ / val, units_ };
                    }

                    // compute a unit addition and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator+(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::op::convert(other.units_, units_));
                        double ntol = std::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                        return { value_ + cval * other.value_, ntol, units_ };
                    }

                    uncertain_measurement simple_add(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::op::convert(other.units_, units_));
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return { value_ + cval * other.value_, ntol, units_ };
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator-(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::op::convert(other.units_, units_));
                        double ntol = std::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                        return { value_ - cval * other.value_, ntol, units_ };
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the simple uncertainty summation method
                    uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                        auto cval = math::op::convert(other.units_, units_);
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return { value_ - cval * other.value_, ntol, units_ };
                    }

                    uncertain_measurement operator+(const measurement& other) const {
                        auto cval = other.value_as(units_);
                        return { value_ + cval, uncertainty_, units_ };
                    }

                    uncertain_measurement operator-(const measurement& other) const {
                        auto cval = other.value_as(units_);
                        return { value_ - cval, uncertainty_, units_ };
                    }

                    // comparison operators 
                    bool operator==(const measurement& other) const {
                        auto val = other.value_as(units_);
                        if (uncertainty_ == 0.0F) { return (value_ == val) ? true : math::op::compare_round_equals(value_, val); }
                        return (val >= (value_ - uncertainty_) && val <= (value_ + uncertainty_));
                    }

                    bool operator>(const measurement& other) const {
                        return value_ > other.value_as(units_);
                    }

                    bool operator<(const measurement& other) const {
                        return value_ < other.value_as(units_);
                    }

                    bool operator>=(const measurement& other) const {
                        auto val = other.value_as(units_);
                        return (value() >= val) ? true : operator==(measurement(val, units_));
                    }

                    bool operator<=(const measurement& other) const {
                        auto val = other.value_as(units_);
                        return (value() <= val) ? true : operator==(measurement(val, units_));
                    }

                    bool operator!=(const measurement& other) const {
                        return !operator==(other);
                    }

                    bool operator==(const uncertain_measurement& other) const {
                        auto zval = simple_subtract(other);
                        return (zval == measurement(0.0, units_));
                    }

                    bool operator>(const uncertain_measurement& other) const {
                        return value_ > other.value_as(units_);
                    }

                    bool operator<(const uncertain_measurement& other) const {
                        return value_ < other.value_as(units_);
                    }

                    bool operator>=(const uncertain_measurement& other) const {
                        auto zval = simple_subtract(other);
                        return (zval.value_ >= 0.0F) ? true :
                                                    (zval == measurement(0.0, units_));
                    }

                    bool operator<=(const uncertain_measurement& other) const {
                        auto zval = simple_subtract(other);
                        return (zval.value_ <= 0.0F) ? true : (zval == measurement(0.0, units_));
                    }

                    bool operator!=(const uncertain_measurement& other) const {
                        return !operator==(other);
                    }


                    // =============================================                                                                                         
                    // get and set methods
                    // =============================================  
                    
                    // set the uncertainty
                    constexpr inline void uncertainty(const double& uncert) { uncertainty_ = uncert; }

                    // set the uncertainty
                    uncertain_measurement& uncertainty(double newUncertainty) {
                        uncertainty_ = newUncertainty;
                        return *this;
                    }

                    // set the uncertainty
                    uncertain_measurement& uncertainty(const measurement& newUncertainty) {
                        uncertainty_ = newUncertainty.value_as(units_);
                        return *this;
                    }

                    // get the uncertainty as a separate measurement
                    measurement uncertainty_measurement() const { return { uncertainty(), units_ }; }

                    // // cast operator to a measurement
                    operator measurement() const { return { value(), units_ }; }

                    // get the numerical value 
                    constexpr double value() const { return value_; }

                    // get the numerical value of the uncertainty
                    constexpr double uncertainty() const { return uncertainty_; }

                    // get the underlying units value
                    units::unit units() const { return units_; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? value_ : math::op::convert(value_, units_, desired_units);
                    }

                    // get the numerical value of the uncertainty as a particular unit
                    double uncertainty_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? uncertainty_ : math::op::convert(uncertainty_, units_, desired_units);
                    }

                    // convert a unit to have a new base
                    uncertain_measurement convert_to(const units::unit& newUnits) const {
                        auto cval = math::op::convert(units_, newUnits);
                        return { cval * value_, uncertainty_ * cval, newUnits };
                    }

                    // print the uncertain measurement
                    inline void print() const { 
                        std::cout << value() << " ± " << uncertainty() << " "; 
                        units_.print(); 
                    }

            }; // class uncertain_measurement

        } // namespace measurements

    } // namespace physics
    

    namespace math {

        // namespace defining some usefull operation
        namespace op {

            inline physics::measurements::measurement operator*(double val, const physics::measurements::measurement& meas) { return meas * val; }

            inline physics::measurements::measurement operator*(double val, const physics::units::unit& unit_base) { return { val, unit_base }; }

            inline physics::measurements::measurement operator*(const physics::units::unit& unit_base, double val) { return { val, unit_base }; }

            inline physics::measurements::fixed_measurement operator*(double v1, const physics::measurements::fixed_measurement& v2) { return {v1 * v2.value(), v2.units()}; }

            inline physics::measurements::fixed_measurement operator*(const physics::measurements::fixed_measurement& v2, double v1) { return {v1 * v2.value(), v2.units()}; }

            inline physics::measurements::uncertain_measurement operator*(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) { return v2.operator*(v1); }

            inline physics::measurements::uncertain_measurement operator*(double v1, const physics::measurements::uncertain_measurement& v2) { return v2.operator*(v1); }

            inline physics::measurements::uncertain_measurement operator*(const physics::measurements::uncertain_measurement& v2, double v1) { return v2.operator*(v1); }


            inline physics::measurements::measurement operator/(double val, const physics::measurements::measurement& meas) {
                return { val / meas.value(), meas.units().inv() };
            }

            inline physics::measurements::measurement operator/(const physics::measurements::measurement& meas, double val) {
                return { val / meas.value(), meas.units().inv() };
            }
            
            inline physics::measurements::measurement operator/(double val, const physics::units::unit& unit_base) {
                return { val, unit_base.inv() };
            }
            
            inline physics::measurements::measurement operator/(const physics::units::unit& unit_base, double val) {
                return { 1.0 / val, unit_base };
            }

            inline physics::measurements::fixed_measurement operator/(double v1, const physics::measurements::fixed_measurement& v2) {
                return {v1 / v2.value(), v2.units().inv()};
            }                
            
            inline physics::measurements::fixed_measurement operator/(const physics::measurements::fixed_measurement& v2, double v1) {
                return {v1 / v2.value(), v2.units().inv()};
            }

            inline physics::measurements::uncertain_measurement operator/(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
                double ntol = v2.uncertainty() / v2.value();
                double nval = v1.value() / v2.value();
                return physics::measurements::uncertain_measurement(nval, nval * ntol, v1.units() / v2.units());
            }

            inline physics::measurements::uncertain_measurement operator/(double v1, const physics::measurements::uncertain_measurement& v2) {
                double ntol = v2.uncertainty() / v2.value();
                double nval = v1 / v2.value();
                return physics::measurements::uncertain_measurement(nval, nval * ntol, v2.units().inv());
            }

            inline physics::measurements::uncertain_measurement operator/(int v1, const physics::measurements::uncertain_measurement& v2) {
                double ntol = v2.uncertainty() / v2.value();
                double nval = static_cast<double>(v1) / v2.value();
                return { nval, nval * ntol, v2.units().inv() };
            }


            inline physics::measurements::fixed_measurement operator+(double v1, const physics::measurements::fixed_measurement& v2) {
                return {v1 + v2.value(), v2.units()};
            }

            inline physics::measurements::fixed_measurement operator-(double v1, const physics::measurements::fixed_measurement& v2) {
                return {v1 - v2.value(), v2.units()};
            }
            
            inline physics::measurements::uncertain_measurement operator+(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
                double cval = math::op::convert(v2.units(), v1.units());
                double ntol = v2.uncertainty() * cval;
                return physics::measurements::uncertain_measurement(v1.value() + cval * v2.value(), ntol, v1.units());
            }

            inline physics::measurements::uncertain_measurement operator-(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
                double cval = math::op::convert(v2.units(), v1.units());
                double ntol = v2.uncertainty() * cval;
                return physics::measurements::uncertain_measurement(v1.value() - cval * v2.value(), ntol, v1.units());
            }


            physics::measurements::measurement pow(const physics::measurements::measurement& meas, int power) {
                return physics::measurements::measurement{ math::op::pow(meas.value(), power), meas.units().pow(power) };
            }

            physics::measurements::fixed_measurement pow(const physics::measurements::fixed_measurement& meas, int power) {
                return physics::measurements::fixed_measurement{math::op::pow(meas.value(), power), meas.units().pow(power)};
            }

            physics::measurements::uncertain_measurement pow(const physics::measurements::uncertain_measurement& meas, int power) {
                auto new_value = math::op::pow(meas.value(), power);
                auto new_tol = ((power >= 0) ? power : -power) * new_value * meas.uncertainty() / meas.value();
                return physics::measurements::uncertain_measurement(new_value, new_tol, meas.units().pow(power));                    
            }

            inline physics::measurements::measurement square(const physics::measurements::measurement& meas) { return pow(meas, 2); }
            
            inline physics::measurements::fixed_measurement square(const physics::measurements::fixed_measurement& meas) { return pow(meas, 2); }

            inline physics::measurements::uncertain_measurement square(const physics::measurements::uncertain_measurement& meas) { return pow(meas, 2); }

            inline physics::measurements::measurement cube(const physics::measurements::measurement& meas) { return pow(meas, 3); }
            
            inline physics::measurements::fixed_measurement cube(const physics::measurements::fixed_measurement& meas) { return pow(meas, 3); }
            
            inline physics::measurements::uncertain_measurement cube(const physics::measurements::uncertain_measurement& meas) { return pow(meas, 3); }


            physics::measurements::measurement root(const physics::measurements::measurement& meas, int power) {
                return physics::measurements::measurement(math::op::root(meas.value(), power), math::op::root(meas.units(), power));
            }

            physics::measurements::fixed_measurement root(const physics::measurements::fixed_measurement& meas, int power) {
                return physics::measurements::fixed_measurement(math::op::root(meas.value(), power), math::op::root(meas.units(), power));
            }

            physics::measurements::uncertain_measurement root(const physics::measurements::uncertain_measurement& um, int power) {
                auto new_value = math::op::root(um.value(), power);
                auto new_tol = new_value * um.uncertainty() / (static_cast<double>((power >= 0) ? power : -power) * um.value());
                return physics::measurements::uncertain_measurement(new_value, new_tol, math::op::root(um.units(), power));
            }

            inline physics::measurements::measurement sqrt(const physics::measurements::measurement& meas) { return root(meas, 2); }
            
            inline physics::measurements::fixed_measurement sqrt(const physics::measurements::fixed_measurement& meas) { return root(meas, 2); }

            inline physics::measurements::uncertain_measurement sqrt(const physics::measurements::uncertain_measurement& meas) { return root(meas, 2); }            

            inline physics::measurements::measurement cbrt(const physics::measurements::measurement& meas) { return root(meas, 3); }

            inline physics::measurements::fixed_measurement cbrt(const physics::measurements::fixed_measurement& meas) { return root(meas, 3); }

            inline physics::measurements::uncertain_measurement cbrt(const physics::measurements::uncertain_measurement& meas) { return root(meas, 3); }


            bool operator==(double val, const physics::measurements::fixed_measurement& v2) {
                return v2 == val;
            }

            bool operator!=(double val, const physics::measurements::fixed_measurement& v2) {
                return v2 != val;
            }

            constexpr bool operator>(double val, const physics::measurements::fixed_measurement& v2) {
                return val > v2.value();
            }

            constexpr bool operator<(double val, const physics::measurements::fixed_measurement& v2) {
                return val < v2.value();
            }

            bool operator>=(double val, const physics::measurements::fixed_measurement& v2) {
                return (val >= v2.value()) ? true : (v2 == val);
            }

            bool operator<=(double val, const physics::measurements::fixed_measurement& v2) {
                return (val <= v2.value()) ? true : (v2 == val);
            }

            bool operator==(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return v2 == other;
            }

            bool operator!=(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return v2 != other;
            }

            constexpr bool operator>(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return other.value() > v2.value();
            }

            constexpr bool operator<(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return other.value() < v2.value();
            }

            bool operator>=(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return (other > v2) ? true : (v2 == other);
            }

            bool operator<=(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return (other < v2) ? true : (v2 == other);
            }

        } // namespace op

    } // namespace math 

    
    namespace physics {

        // namespace for keeping track of an object in a 3D system
        namespace position {

            class coordinate {

                protected: 

                    // =============================================
                    // class members
                    // =============================================
                
                    measurements::fixed_measurement m_coordinate; 
                    

                public:  

                    // =============================================
                    // constructors and destructor
                    // =============================================

                    coordinate(const double& coord, const units::unit& unit) : m_coordinate(coord, unit) {}
                    
                    coordinate(const units::unit& unit, const double& coord) : m_coordinate(coord, unit) {}

                    coordinate(const measurements::fixed_measurement& coord) : m_coordinate(coord) {}

                    coordinate() = default; 


                    // =============================================
                    // set, get and print methods
                    // =============================================
                    
                    constexpr void value(const double& x) { m_coordinate.value(x); }

                    double value() const { return m_coordinate.value(); }

                    measurements::fixed_measurement get_coordinate_measurement() const { return m_coordinate; }

                    void print() const { m_coordinate.print(); }

            }; // class coordinate


            class position {

                private: 
                    
                    // =============================================
                    // class members
                    // =============================================
                    
                    std::vector<coordinate> m_position; 
                

                public: 

                    // =============================================
                    // constructors
                    // =============================================

                    position(const coordinate& x, const coordinate& y, const coordinate& z) {
                        m_position.push_back(x); 
                        m_position.push_back(y);
                        m_position.push_back(z); 
                    }

                    position(const measurements::fixed_measurement& x, const measurements::fixed_measurement& y, const measurements::fixed_measurement& z) {
                        m_position.push_back(x); 
                        m_position.push_back(y);
                        m_position.push_back(z);
                    }

                    position(const std::vector<coordinate>& pos) : m_position{pos} {}

                    position(const position& pos) : position(pos.get_coordinates()) {}


                    // =============================================
                    // set, get and print methods
                    // =============================================

                    void set_position(const std::vector<coordinate>& pos) { m_position = pos; }

                    void coordinate_x(const coordinate& x) { m_position[0] = x; }
                    
                    void coordinate_y(const coordinate& y) { m_position[1] = y; }
                    
                    void coordinate_z(const coordinate& z) { m_position[2] = z; }

                    void value_x(const double& x) { m_position[0].value(x); }

                    void value_y(const double& y) { m_position[1].value(y); }

                    void value_z(const double& z) { m_position[2].value(z); }

                    std::vector<coordinate> get_coordinates() const { return m_position; }

                    coordinate coordinate_x() const { return m_position[0]; }

                    coordinate coordinate_y() const { return m_position[1]; }

                    coordinate coordinate_z() const { return m_position[2]; }

                    double value_x() const { return m_position[0].value(); }

                    double value_y() const { return m_position[1].value(); }

                    double value_z() const { return m_position[2].value(); }

                    measurements::fixed_measurement magnitude() const { 
                        return math::op::sqrt(math::op::square(m_position[0].get_coordinate_measurement()) +
                                                        math::op::square(m_position[1].get_coordinate_measurement()) +
                                                        math::op::square(m_position[2].get_coordinate_measurement()));
                    }

                    measurements::fixed_measurement distance(const std::vector<coordinate>& pos) const {        
                        return math::op::sqrt(math::op::square(pos[0].get_coordinate_measurement() - m_position[0].get_coordinate_measurement()) +
                                                        math::op::square(pos[1].get_coordinate_measurement() - m_position[1].get_coordinate_measurement()) +
                                                        math::op::square(pos[2].get_coordinate_measurement() - m_position[2].get_coordinate_measurement()));
                    }                    

                    measurements::fixed_measurement distance(const position& pos) const {        
                        return math::op::sqrt(math::op::square(pos.coordinate_x().get_coordinate_measurement() - m_position[0].get_coordinate_measurement()) +
                                                    math::op::square(pos.coordinate_y().get_coordinate_measurement() - m_position[1].get_coordinate_measurement()) +
                                                    math::op::square(pos.coordinate_z().get_coordinate_measurement() - m_position[2].get_coordinate_measurement()));
                    }       
                    
                    double rho() { return math::op::sqrt(math::op::square(value_x()) + math::op::square(value_y())); }

                    double phi() const { return atan2(value_y(), value_x()); }     

                    double phi(const std::vector<coordinate>& pos) const { 
                        return atan2(pos[1].value() - value_y(), pos[0].value() - value_x()); 
                    }

                    double phi(const position& pos) const { 
                        return atan2(pos.value_y() - value_y(), pos.value_x() - value_x()); 
                    }

                    double theta() const { return acos(value_z() / magnitude().value()); }
            
                    double theta(const std::vector<coordinate>& pos) const { 
                        return acos((pos[2].value() - value_z()) / distance(pos).value()); 
                    }

                    double theta(const position& pos) const { 
                        return acos((pos.value_z() - value_z()) / distance(pos).value()); 
                    }
    
                    // std::vector<coordinate> direction() const {
                    //     return {cos(phi()), sin(phi()), m_position[2].coordinate() / magnitude()}};
                    // } 

                    // std::vector<coordinate> direction(const std::vector<coordinate>& pos1) const {
                    //     return {cos(phi(pos1)), sin(phi(pos1)), (pos1[2].coordinate() - m_position[2].coordinate()) / distance(pos1)};
                    // } 

                    // std::vector<coordinate> direction(const position& pos1) const {
                    //     return {cos(phi(pos1)), sin(phi(pos1)), (pos1.coordinate_z().coordinate() - m_position[2].coordinate()) / distance(pos1)};
                    // } 

                    void print() const {
                        std::cout << "- position = { ";
                        for (auto i : m_position) { 
                            std::cout << "\t["; 
                            i.coordinate::print(); 
                            std::cout << "]";
                        }
                        std::cout << "  }" << std::endl; 
                    }


            }; // class position


            class velocity {

                private: 
                    
                    // =============================================
                    // class members
                    // =============================================
                    
                    std::vector<coordinate> m_velocity; 
                

                public:  

                    // =============================================
                    // constructors
                    // =============================================

                    velocity(const coordinate& x, const coordinate& y, const coordinate& z) {
                        m_velocity.push_back(x); 
                        m_velocity.push_back(y);
                        m_velocity.push_back(z); 
                    }

                    velocity(const measurements::fixed_measurement& x, const measurements::fixed_measurement& y, const measurements::fixed_measurement& z) {
                        m_velocity.push_back(x); 
                        m_velocity.push_back(y);
                        m_velocity.push_back(z);
                    }

                    velocity(const std::vector<coordinate>& vel) : m_velocity{vel} {}

                    velocity(const velocity& vel) : m_velocity{vel.get_coordinates()} {}


                    // =============================================
                    // set, get and print methods
                    // =============================================
                    void set_velocity(const std::vector<coordinate>& pos) { m_velocity = pos; }

                    void coordinate_x(const coordinate& x) { m_velocity[0] = x; }
                    
                    void coordinate_y(const coordinate& y) { m_velocity[1] = y; }
                    
                    void coordinate_z(const coordinate& z) { m_velocity[2] = z; }

                    void value_x(const double& x) { m_velocity[0].value(x); }

                    void value_y(const double& y) { m_velocity[1].value(y); }

                    void value_z(const double& z) { m_velocity[2].value(z); }

                    std::vector<coordinate> get_coordinates() const { return m_velocity; }

                    coordinate coordinate_x() const { return m_velocity[0]; }

                    coordinate coordinate_y() const { return m_velocity[1]; }

                    coordinate coordinate_z() const { return m_velocity[2]; }

                    double value_x() const { return m_velocity[0].value(); }

                    double value_y() const { return m_velocity[1].value(); }

                    double value_z() const { return m_velocity[2].value(); }

                    measurements::fixed_measurement magnitude() const { 
                        return math::op::sqrt(math::op::square(m_velocity[0].get_coordinate_measurement()) +
                                                        math::op::square(m_velocity[1].get_coordinate_measurement()) +
                                                        math::op::square(m_velocity[2].get_coordinate_measurement()));
                    }
                
                    double phi() const { return atan2(value_y(), value_x()); }     
                    
                    double theta() const { return acos(value_z() / magnitude().value()); }

                    // std::vector<coordinate> direction() const {
                    //     return {cos(phi()), sin(phi()), cos(theta())};
                    // } 

                    void print() const {
                        std::cout << "- velocity = { ";
                            for (auto i : m_velocity) { 
                            std::cout << "\t["; 
                            i.coordinate::print(); 
                            std::cout << "]";
                        }
                        std::cout << "  }" << std::endl; 
                    }

            }; // class velocity
            
        } // namespace position
    
    } // namespace physics


    namespace math {

        namespace functions {

            class function_base {

                public: 
                
                    // =============================================
                    // constructor and destructor
                    // =============================================     
                
                    function_base() {}

                    ~function_base() {}

                
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
                
                    virtual void print() const = 0;
                
                    
                    // =============================================
                    // extra methods
                    // =============================================
                
                    int signum(const double& x) const { return (x == 0. ? 0. : (x > 0 ? 1. : -1)); }
                
            }; // class function_base


            class functor : public function_base {

                private:

                    // =============================================
                    // class members
                    // =============================================
                
                    char m_op; 
                    function_base * m_f; 
                    function_base * m_g;
                
                
                public:

                    // =============================================
                    // constructor and destructor
                    // =============================================     
                
                    functor(const char& op, function_base& f, function_base& g) {
                        m_f = &f; 
                        m_g = &g;
                        m_op = op;
                    }
                    
                    ~functor() {}
                
                
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
                                return std::pow(m_f->eval(x), m_g->eval(x));

                            case 'c':
                                return m_f->eval(m_g->eval(x));
                        }
                        return NAN;
                    }

                    
                    // =============================================
                    // print methods
                    // =============================================
                
                    void print() const override {}
                    
            };


            class quadratic : public function_base {

                private: 

                    // =============================================
                    // class members
                    // =============================================
                    
                    double m_a, m_b, m_c, m_delta; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================     
            
                    explicit quadratic(const double& a, const double& b, const double& c) noexcept : 
                        m_a{a}, m_b{b}, m_c{c}, m_delta{pow(m_b, 2) - 4 * m_a * m_c} {}
                    
                    ~quadratic() = default; 

            
                    // =============================================
                    // set methods
                    // =============================================
            
                    constexpr inline void a(const double& a) { m_a = a; } 

                    constexpr inline void b(const double& b) { m_b = b; }

                    constexpr inline void c(const double& c) { m_c = c; }  
            
            
                    // =============================================
                    // get methods
                    // =============================================

                    constexpr inline double a() const { return m_a; }

                    constexpr inline double b() const { return m_b; } 

                    constexpr inline double c() const { return m_c; } 

                    constexpr inline double delta() const { return m_delta; } 

                    std::vector<double> roots() const {
                        if (m_delta == 0) return { - m_b / (2 * m_a), - m_b / (2 * m_a) };
                        else if (m_delta > 0) return { (- m_b - std::sqrt(m_delta)) / (2 * m_a), (- m_b + std::sqrt(m_delta)) / (2 * m_a) };
                        else std::cout << std::endl << "There are not real solutions..." << std::endl << std::endl; // note: create a complex class then come back here
                        exit(-11);  
                    }


                    // =============================================
                    // eval methods
                    // =============================================

                    constexpr inline double eval_Horner(const double& x) const { return m_c + x * (m_b + x * m_a); }

                    constexpr inline double eval(const double& x) const override { return m_a * math::op::square(x) + m_b * x + m_c; }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print() const override {
                        std::cout << " y = "; 
                        if (m_a != 1) std::cout << m_a; 
                        std::cout << "x^2 "; 
                        if (m_b != 0) {
                            if (m_b != 1) std::cout << "+ " << m_b << "x "; 
                            else std::cout << "+ x ";                       
                        }
                        if (m_c != 0) { std::cout << "+ " << m_c << std::endl; }
                        else std::cout << std::endl;
                    }                
            
            }; 


            class cubic : public function_base {
                
                private: 

                    // =============================================
                    // class members
                    // =============================================
                    
                    double m_a, m_b, m_c, m_d; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
                
                    cubic(double a, double b, double c, double d) : m_a{a}, m_b{b}, m_c{c}, m_d{d} {}
                
                    ~cubic() {} 
                
                
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

                    void print() const override {
                        std::cout << " y = "; 
                        if (m_a != 1) std::cout << m_a; 
                        std::cout << "x^3 "; 
                        if (m_b != 0) {
                            if (m_b != 1) std::cout << "+ " << m_b << "x^2 "; 
                            else std::cout << "+ x^2 "; 
                        }
                        if (m_c != 0) {
                            if (m_c != 1) std::cout << "+ " << m_c << "x"; 
                            else std::cout << "+ x"; 
                        }
                        if (m_d != 0) std::cout << "+ " << m_d << std::endl;
                        else std::cout << std::endl;  
                    }
            
            }; 


            class squareRoot : public function_base {
                
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
                
                    squareRoot(double c1 = 1) : m_c1{c1} {}
                
                    ~squareRoot() {} 
                
                
                    // =============================================
                    // set & get methods
                    // =============================================
                
                    void set_c1(double c1) { m_c1 = c1; }

                    double get_c1() const { return m_c1; }

                
                    // =============================================
                    // eval methods
                    // =============================================

                    double eval(const double& x) const override { return m_c1 * pow(x, 0.5); }
                
                    double eval(const function_base& f, double x) const { return m_c1 * pow(f.eval(x), 0.5); }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1; 
                        std::cout << "x^(1/2)" << std::endl;
                    }        
                
                
            };


            class cubicRoot : public function_base {
                
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
                
                    cubicRoot(double c1 = 1) : m_c1{c1} {}
                
                    ~cubicRoot() {} 
                
                
                    // =============================================
                    // set & get methods
                    // =============================================
                
                    void set_c1(double c1) { m_c1 = c1; }

                    double get_c1() const { return m_c1; }

                
                    // =============================================
                    // eval methods
                    // =============================================

                    double eval(const double& x) const override { return m_c1 * pow(x, 1. / 3.); }
                
                    double eval(const function_base& f, double x) const { return m_c1 * pow(f.eval(x), 1. / 3.); }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1; 
                        std::cout << "x^(1/3)" << std::endl;
                    }        
                
            };


            class exponential : public function_base {
            
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
            
                    exponential(double base = physics::constants::e, double c1 = 1, double c2 = 1) : m_base{base}, m_c1{c1}, m_c2{c2} {}
                
                    ~exponential() {} 
            

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
            
                    double eval(const function_base& f, double x) const { return m_c1 * pow(m_base, m_c2 * f.eval(x)); }
                
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        if (m_base != physics::constants::e) std::cout << m_base << "^"; 
                        else std::cout << "e^(";
                        if (m_c2 != 1) std::cout << m_c2; 
                        std::cout << "x)" << std::endl; 
                    }

            }; 


            class logarithm : public function_base {
            
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
            
                    logarithm(double base = physics::constants::e, double c1 = 1, double c2 = 1) : m_base{base}, m_c1{c1}, m_c2{c2} {}
                
                    ~logarithm() {} 
            

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
            
                    double eval(const double& x) const override { return m_c1 * std::log(m_c2 * x) / std::log(m_base); }
                    
                    double eval(const function_base& f, double x) const { return m_c1 * std::log(m_c2 * f.eval(x)) / std::log(m_base); }
            
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "log";
                        if (m_base != physics::constants::e) std::cout << "_( " << m_base << ")";
                        if (m_c2 != 1) std::cout << "(" << m_c2 << "x)" << std::endl; 
                        else std::cout << "(x)" << std::endl;
                    }                

            }; 


            class sine : public function_base {
                
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
            
                    sine(double c1 = 1, double c2 = 1) : m_c1{c1}, m_c2{c2} {} 
                
                    ~sine() {} 

                
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
            
                    double eval(const double& x) const override { return m_c1 * std::sin(m_c2 * x); }
            
                    double eval(const function_base& f, double x) const { return m_c1 * std::sin(m_c2 * f.eval(x)); }
            
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "sin(";
                        if (m_c2 != 1) std::cout << m_c2; 
                        std::cout << "x)" << std::endl;
                    }

            }; 


            class cosine : public function_base {
                
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
            
            
                    cosine(double c1 = 1, double c2 = 1) : m_c1{c1}, m_c2{c2} {} 
                
                    ~cosine() {} 

                
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
            
                    double eval(const double& x) const override { return m_c1 * std::cos(m_c2 * x); }
            
                    double eval(const function_base& f, double x) const { return m_c1 * std::cos(m_c2 * f.eval(x)); }
            
            
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "cos(";
                        if (m_c2 != 1) std::cout << m_c2; 
                        std::cout << "x)" << std::endl;
                    }

            }; 


            class tangent : public function_base {

                private: 
                    // =============================================
                    // class members
                    // =============================================   
                    
                    // y = c1 * tan (c2 * x)
                    double m_c1, m_c2; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
            
                    tangent(double c1 = 1, double c2 = 1) : m_c1{c1}, m_c2{c2} {} 
                
                    ~tangent() {} 

                
                    // =============================================
                    // set methods
                    // =============================================
            
                    void set_c1(double c1) { m_c1 = c1; }
                
                    void set_c2(double c2) { m_c2 = c2; }
                
                
                    // =============================================
                    // get methods
                    // =============================================

                    double inline get_c1() const { return m_c1; }

                    double inline get_c2() const { return m_c2; }
                

                    // =============================================
                    // eval methods
                    // =============================================
            
                    double eval(const double& x) const override { return m_c1 * std::tan(m_c2 * x); }
            
                    double eval(const function_base& f, double x) const { return m_c1 * std::tan(m_c2 * f.eval(x)); }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "tan(";
                        if (m_c2 != 1) std::cout << m_c2; 
                        std::cout << "x)" << std::endl;
                    }

            }; // class tangent

        } // namespace functions
        

        namespace statistics {

            // mean
            template <typename T, typename IT>
            T mean(IT begin, IT end) {
                if (end - begin == 0) return 0; 
                return std::accumulate(begin, end, 0.0) / (end - begin); 
            }

            template <typename T, typename K>
            T mean(const std::vector<K>& v) {
                T accu{}; 
                if (v.size() == 0) return accu;
                for (T x : v) accu += x;
                return accu / v.size(); 
            }

            // median
            template <typename T, typename IT>
            T median(IT begin, IT end) {
                sort(begin, end); 
                if ((end - begin) % 2 != 0) return *(begin + (end - begin) / 2);
                else return ((*(begin + (end - begin) / 2) + *(begin + (end - begin) / 2 - 1)) / 2);
            }

            template <typename T, typename K>
            T median(std::vector<K> v) {
                sort(v.begin(), v.end()); 
                if (v.size()%2 != 0) return v[v.size() / 2];
                else return (v[v.size() / 2] + v[(v.size() / 2)-1]) / 2; 
            }

            // variance
            template <typename T, typename IT>
            T var(IT begin, IT end) {
                if (end - begin == 0) return 0; 
                T accu{}; 
                for (IT i{begin}; i < end; i++) accu += pow(*i - mean<T>(begin, end), 2);
                return accu / (end - begin);
            }

            template <typename T, typename K>
            T var(std::vector<K> v) {
                T accu{}; 
                if (v.size() == 0) return accu;
                for (auto x : v) accu += pow(x - mean<T>(v), 2);
                return accu / v.size();
            }

            // standard deviation
            template <typename T, typename IT>
            T sd(IT begin, IT end) {
                return sqrt(var<T>(begin, end));
            }

            template <typename T, typename K>
            T sd(std::vector<K> v) {
                return sqrt(var<T>(v));
            }

            // standard deviation of mean
            template <typename T, typename IT>
            T sdom(IT begin, IT end) {
                return sd(begin, end) / sqrt(end - begin);
            }

            template <typename T, typename K>
            T sdom(std::vector<K> v) {
                return sd(v) / sqrt(v.size());
            }

            // chi squared
            template <typename T, typename IT, typename J>
            T chi_sq(IT begin, IT end, J expected_value) {
                T accu{}; 
                for (auto x : end) accu += pow(x - expected_value, 2) / sd(begin, end); 
                return accu; 
            }

            template <typename T, typename K, typename J>
            T chi_sq(std::vector<K> v, J expected_value) {
                T accu{}; 
                for (auto x : v) accu += pow(x - expected_value, 2) / sd(v); 
                return accu; 
            }
            
        } // namespace statistics


        class random_generator {

            private: 

                // =============================================
                // class members
                // =============================================        
                
                unsigned int m_a, m_c, m_m, m_seed;


            public: 
        
                // =============================================
                // constructors
                // =============================================
                
                random_generator() { set_a(1664525); set_c(1013904223); set_m(pow(2, 31)); }

                random_generator(unsigned int seed) : random_generator() { m_seed = seed; }


                // =============================================
                // set and get methods
                // =============================================

                void set_a(unsigned int a) { m_a = a; }
                
                void set_c(unsigned int c) { m_c = c; }
                
                void set_m(unsigned int m) { m_m = m; }
                
                void set_seed(unsigned int seed) { m_seed = seed; }

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
                    return - std::log(1 - rand()) / mean; 
                }

                double gauss_box_muller(double mean, double sigma) {
                    double s{rand()}, t{rand()}, x{};
                    x = std::sqrt(-2 * log(s)) * cos(2 * M_PI * t);
                    return mean + sigma * x;
                }

                double gauss_accept_reject(double mean, double sigma) {
                    double x{}, y{}, g{}; 
                    while (true) {
                        x = rand(-5., 5.); 
                        y = rand(); 
                        g = exp(pow(x, 2) / 2); 
                        if (y <= g) break;
                    }
                    return mean + x * sigma;
                }

        }; // class random_generator


        class integral {

            protected: 

                // =============================================
                // class members
                // =============================================

                double m_sum, m_integral, m_old_integral, m_error;  

                double m_a, m_b, m_h;

                unsigned int m_steps;

                int m_sign; 

                random_generator m_rg;


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

                integral() : m_a{}, m_b{}, m_h{}, m_sign{}, m_sum{},
                    m_integral{}, m_old_integral{}, m_error{} {}
                
                ~integral() = default; 

            
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

                void print_integral(const double& precision = 1.e-6) {
                    std::cout << "\nIntegral of f(x) in [" << m_a << ", " << m_b << "] = " << std::setprecision((int)-log10(precision)) << m_integral << std::endl;
                }

                void print_error() {
                    std::cout << "error = " << m_error << std::endl;
                }        
                
                
                // =============================================
                // integration methods
                // =============================================

                void midpoint(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    for (unsigned int i{}; i < m_steps; i++) { m_sum += (f.eval(m_a + (i + 0.5) * m_h)); }
                    m_integral = m_sign * m_sum * m_h; 
                }

                void midpoint_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-5) {
                    begin_integration(a, b, 1); 
                    double oldintegral;  
                    while (true) {
                        oldintegral = m_integral; 
                        midpoint(m_a, m_b, f, m_steps * 2);
                        error_integral(fabs(m_integral - oldintegral));
                        if (m_error < prec) break;
                    }    
                }
                
                void trapexoid(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 2.);
                    for (unsigned int i{1}; i < m_steps; i++) { m_sum += f.eval(m_a + i * m_h); }
                    m_integral = m_sign * m_sum * m_h; 
                } 

                void trapexoid_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-5) {
                    begin_integration(a, b, 2, f.eval(a) + f.eval(b) / 2. + f.eval((a + b) / 2.)); 
                    double oldintegral;  
                    while (true) {
                        oldintegral = m_integral; 
                        trapexoid(m_a, m_b, f, m_steps * 2);
                        error_integral(fabs(m_integral - oldintegral));
                        if (m_error < prec) break;
                    }
                }

                void simpson(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    if (n % 2 == 0) begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 3.);
                    else begin_integration(a, b, n + 1);  
                    for (unsigned int i{1}; i < m_steps; i++) {
                        m_sum += 2 * (1 + i % 2) * (f.eval(m_a + i * m_h)) / 3.; 
                    }
                    m_integral = m_sign * m_sum * m_h; 
                }

                void simpson_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-5) {
                    begin_integration(a, b, 2, (f.eval(a) + f.eval(b)) / 3.); 
                    double oldintegral; 
                    while (true) {
                        oldintegral = m_integral; 
                        simpson(m_a, m_b, f, m_steps * 2);
                        error_integral(fabs(m_integral - oldintegral));
                        if (m_error < prec) break; 
                    }
                }

                void mean(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    for (unsigned int i{}; i < n; i ++) {
                        m_sum += f.eval(m_rg.rand(a, b)); 
                    }
                    m_integral = (b - a) * m_sum / n; 
                }

                void mean_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-5) {
                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        mean(a, b, f);
                        k.push_back(m_integral); 
                    }
                    double k_mean = sqrt(100) * statistics::sd<double>(k); 
                    unsigned int N = (unsigned int) pow(k_mean / prec, 2); 
                    mean(a, b, f, N); 
                }
        
                void hit_or_miss(const double& a, const double& b, const functions::function_base& f, const double& fmax, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    double x{}, y{}; 
                    unsigned int hits{};
                    for (unsigned int i{}; i < n; i ++) {
                        x = m_rg.rand(a, b); 
                        y = m_rg.rand(0., fmax);  
                        if (y <= f.eval(x)) hits++; 
                    }
                    m_integral = (b - a) * fmax * hits / n; 
                }

                void hit_or_miss_fixed(const double& a, const double& b, const functions::function_base& f, const double& fmax, const double& prec = 1.e-5) {
                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        hit_or_miss(a, b, f, fmax);
                        k.push_back(m_integral); 
                    }
                    double k_mean = sqrt(100) * statistics::sd<double>(k); 
                    unsigned int N = (unsigned int) pow(k_mean / prec, 2); 
                    hit_or_miss(a, b, f, fmax, N); 
                }


        }; // class integral

    } // namespace math


    namespace physics {
        namespace mass {

            class mass : public position::position {

                protected: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::fixed_measurement m_mass; 
                    

                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    mass(const double& mass, const units::unit& mass_unit, const std::vector<position::coordinate>& pos) : 
                        m_mass{mass, mass_unit}, position::position(pos) {}

                    mass(const double& mass, const units::unit& mass_unit, const position::position& pos) : 
                        m_mass{mass, mass_unit}, position::position(pos) {}

                    mass(const double& mass, const std::vector<position::coordinate>& pos) : 
                        m_mass{mass, units::defined::kg}, position::position(pos) {}

                    mass(const double& mass, const position::position& pos) : 
                        m_mass{mass, units::defined::kg}, position::position(pos) {}

                    ~mass() {}


                    // =============================================
                    // set, get and print methods
                    // =============================================

                    void set_mass(const double& mass) { m_mass.value(mass); }
                    
                    double inline get_mass() const { return m_mass.value(); }

                    measurements::fixed_measurement get_mass_measurement() const { return m_mass; }

                    void inline print() const { std::cout << "- mass = " << get_mass() << " udm" << std::endl; }
                
            }; // class mass

        } // namespace mass
        

        namespace charge {
        
            class charge : public position::position {

                protected: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::fixed_measurement m_charge; 
                    

                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    charge(const double& charge, const units::unit& charge_unit, const std::vector<position::coordinate>& pos) : 
                        m_charge{charge, charge_unit}, position::position(pos) {}

                    charge(const double& charge, const units::unit& charge_unit, const position::position& pos) : 
                        m_charge{charge, charge_unit}, position::position(pos) {}

                    charge(const double& charge, const std::vector<position::coordinate>& pos) : 
                        m_charge{charge, units::defined::C}, position::position(pos) {}

                    charge(const double& charge, const position::position& pos) : 
                        m_charge{charge, units::defined::C}, position::position(pos) {}

                    ~charge() {}


                    // =============================================
                    // set, get and print methods
                    // =============================================

                    void set_charge(const double& charge) { m_charge.value(charge); }
                    
                    double inline get_charge() const { return m_charge.value(); }

                    measurements::fixed_measurement get_charge_measurement() const { return m_charge; }

                    void inline print() const { std::cout << "- charge = " << get_charge() << " udm" << std::endl; }
                
            }; // class charge

        } // namespace charge


        namespace time {

            class time {

                public:     

                    // =============================================
                    // class members
                    // =============================================     

                    measurements::fixed_measurement m_time; 


                    // =============================================
                    // constructor and destructor
                    // =============================================   

                    time(const double& time, const units::unit& unit) : m_time(time, unit) {}
                    
                    time(const units::unit& unit, const double& time) : m_time(time, unit) {}

                    time(const measurements::fixed_measurement& time) : m_time(time) {}

                    time(const units::unit& unit) : m_time(0.0, unit) {}

                    time(const double& time) : m_time(time, units::defined::s) {}

                    time() : m_time(0.0, units::defined::s) {}

                    ~time() = default; 

                    
                    // =============================================
                    // time methods
                    // =============================================   

                    double get_time() const { return m_time.value(); }
                                        
                    measurements::fixed_measurement get_time_measurement() const { return m_time; }

                    void increase_time(const double& h) { m_time.value(m_time.value() + h); } 

                    void reset_time() { m_time.value(0); }

                    void print() const { m_time.print(); }

            }; // class time


            class timer : public time {

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

                    timer(const units::unit& unit) : time(0.0, unit) {}

                    timer() : time(0.0, units::defined::s) {}

                    ~timer() = default;

                    
                    // =============================================
                    // timer methods
                    // =============================================   

                    void start() { m_start = std::chrono::system_clock::now(); }

                    void pause() { m_pause = std::chrono::system_clock::now(); }

                    void end() { m_end = std::chrono::system_clock::now(); }

                    double elapsed_time() { 
                        m_elapsed_seconds = m_pause - m_start;
                        return m_elapsed_seconds.count(); 
                    }
                    
                    void print_elapsed_time() { std::cout << "Elapsed time = " << elapsed_time() << " udm" << std::endl; }

            }; // class timer

        } // namespace time

    } // namespace physics

} // namespace physim


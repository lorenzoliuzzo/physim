
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools for computational physics. 
// last updated:    24/07/2022

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional> 
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>


namespace physics {

    namespace tools {

        // namespace defining some usefull constants
        namespace constants {

            constexpr int32_t max_neg(uint32_t n_bits) { return -(int32_t(1U << (n_bits - 1))); }
            constexpr double infinity = std::numeric_limits<double>::infinity();
            constexpr double invalid_conversion = std::numeric_limits<double>::signaling_NaN();

            constexpr double pi = 3.14159265358979323846;

        } // namespace constants


        // verify if the number is NaN 
        static_assert(std::numeric_limits<double>::has_signaling_NaN, "nan is used to signify invalid values");
            
        // verify if the number is infinity
        static_assert(std::numeric_limits<double>::has_infinity, "nan is used to signify invalid values");


        // namespace defining some usefull operation
        namespace op {

            // operator to generate the square power of a value
            template<typename T> constexpr T sqr_power(T a) { return a * a; }

            // operator to generate the cubic power of a value
            template<typename T> constexpr T cbc_power(T a) { return a * a * a; }

            // operator to generate small integer powers of a value(1, 0, -1)
            template<typename T> constexpr T power_const_small(T val, int power) { 
                return (power == 1) ? val : ((power == -1) ? T(1.0) / val : T(1.0));
            }

            // operator to generate an integer power of a number
            template<typename T> constexpr T power_const(T val, int power) {
                return (power > 1) ? sqr_power(power_const(val, power / 2)) * (power % 2 == 0 ? T(1.0) : val) :
                    (power < -1) ? T(1.0) / (sqr_power(power_const(val, (-power) / 2)) * ((-power) % 2 == 0 ? T(1.0) : val)) :
                    power_const_small(val, power);
            }

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

        }

        
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

                public:

                    // the seven SI base units 
                    enum base {
                        Meter = 0,
                        Second = 1,
                        Kilogram = 2,
                        Ampere = 3,
                        Kelvin = 4,
                        Mole = 5,
                        Candela = 6,  
                    };

                    static constexpr uint32_t bits[14] = { bitwidth::meter,
                                                           bitwidth::second,
                                                           bitwidth::kilogram,
                                                           bitwidth::ampere,
                                                           bitwidth::kelvin,
                                                           bitwidth::mole,
                                                           bitwidth::candela
                                                        };  


                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // constructor from powers
                    constexpr unit_data(int meters, 
                                        int seconds, 
                                        int kilograms,
                                        int amperes, 
                                        int kelvins, 
                                        int moles, 
                                        int candelas) :
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
                    constexpr unit_data pow(int power) const { 
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
                    constexpr unit_data root(int power) const {
                        return (has_valid_root(power)) ? unit_data( meter_ / power,
                                                                    second_ / power,
                                                                    kilogram_ / power,
                                                                    ampere_ / power,
                                                                    kelvin_ / power,
                                                                    mole_ / power,
                                                                    candela_ / power) : 
                                                                        unit_data(nullptr);
                    }
                    
                    
                    // =============================================
                    // check methods
                    // ============================================= 

                    // comparison operators
                    constexpr bool operator==(const unit_data& other) const {
                        return equivalent_non_counting(other) && mole_ == other.mole_;         
                    }

                    constexpr bool operator!=(const unit_data& other) const {
                        return !(*this == other);
                    }

                    // check if the units have the same base unit 
                    constexpr bool has_same_base(const unit_data& other) const {
                        return equivalent_non_counting(other) && mole_ == other.mole_;
                    }

                    // check equivalence for non-counting base units
                    constexpr bool equivalent_non_counting(const unit_data& other) const {
                        return meter_ == other.meter_ && second_ == other.second_ &&
                            kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                            candela_ == other.candela_ && kelvin_ == other.kelvin_;
                    }
                    
                    // check if the unit is empty
                    constexpr bool empty() const {
                        return meter_ == 0 && second_ == 0 && kilogram_ == 0 &&
                            ampere_ == 0 && candela_ == 0 && kelvin_ == 0 && mole_ == 0;
                    }
                    

                    // =============================================
                    // get powers
                    // =============================================
                    
                    constexpr int meter() const { return meter_; }

                    constexpr int second() const { return second_; }

                    constexpr int kg() const { return kilogram_; }
                                        
                    constexpr int ampere() const { return ampere_; }
                    
                    constexpr int kelvin() const { return kelvin_; }
                    
                    constexpr int mole() const { return mole_; }
                    
                    constexpr int candela() const { return candela_; }

                    // Get the number of different base units used
                    constexpr int unit_type_count() const {
                        return ((meter_ != 0) ? 1 : 0) + ((second_ != 0) ? 1 : 0) +
                            ((kilogram_ != 0) ? 1 : 0) + ((ampere_ != 0) ? 1 : 0) +
                            ((candela_ != 0) ? 1 : 0) + ((kelvin_ != 0) ? 1 : 0) +
                            ((mole_ != 0) ? 1 : 0);
                    }


                private: 

                    signed int meter_ : bitwidth::meter;
                    signed int second_ : bitwidth::second;  
                    signed int kilogram_ : bitwidth::kilogram;
                    signed int ampere_ : bitwidth::ampere;
                    signed int candela_ : bitwidth::candela;  
                    signed int kelvin_ : bitwidth::kelvin;
                    signed int mole_ : bitwidth::mole;

                    // check if the base_unit has a valid root
                    constexpr bool has_valid_root(int power) const {
                        return meter_ % power == 0 && second_ % power == 0 &&
                            kilogram_ % power == 0 && ampere_ % power == 0 &&
                            candela_ % power == 0 && kelvin_ % power == 0 &&
                            mole_ % power == 0;
                    }      
                    
                    // to handle a few weird operations that operate on square_root Hz
                    constexpr int root_Hertz_modifier(int power) const {
                        return (second_ * power == 0 || power % 2 != 0) ? 0 :
                            (power / 2) * ((second_ < 0) || (power < 0) ? 9 : -9);
                    }    

            }; // class unit_base


            // class defining a basic unit module with double precision on the multiplier
            class unit {
                
                private:

                    // =============================================
                    // class members
                    // ============================================= 
                    
                    unit_data base_units_{0, 0, 0, 0, 0, 0, 0};

                    double multiplier_{1.0};  


                public:
                    
                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // default constructor
                    constexpr unit() noexcept {};

                    // constructor from base_unit 
                    explicit constexpr unit(const unit_data& base_unit) noexcept : 
                        base_units_(base_unit) {};

                    // constructor from unit with a double 
                    constexpr unit(const unit& other, double mult) noexcept :
                        unit(other.base_units(), mult * other.multiplier()) {};

                    // constructor from double with an unit
                    constexpr unit(double mult, const unit& other) noexcept :
                        unit(other.base_units(), mult * other.multiplier()) {};
                    
                    // constructor from base_unit with a double
                    constexpr unit(const unit_data& base_unit, double mult) noexcept :
                        base_units_(base_unit), multiplier_(mult) {}

                    // constructor from double with a base_unit
                    constexpr unit(double mult, const unit_data& base_unit) noexcept :
                        base_units_(base_unit), multiplier_(mult) {}
                        

                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a unit
                    constexpr unit inv() const { return { base_units_.inv(), 1.0 / multiplier() }; }

                    // multiply with another unit
                    constexpr unit operator*(const unit& other) const {
                        return { base_units_ * other.base_units_, multiplier() * other.multiplier() };
                    }

                    // division operator
                    constexpr unit operator/(const unit& other) const {
                        return { base_units_ / other.base_units_, multiplier() / other.multiplier() };
                    }

                    // take a unit to a power
                    constexpr unit pow(int power) const {
                        return { base_units_.pow(power), op::power_const(multiplier_, power) };
                    }


                    // =============================================
                    // checks
                    // ============================================= 

                    // equality operator
                    bool operator==(const unit& other) const {
                        if (base_units_ != other.base_units_ ) { return false; }
                        if (multiplier_ == other.multiplier_) { return true; }
                        return op::compare_round_equals(multiplier_, other.multiplier_);
                    }

                    // equality operator
                    bool operator!=(const unit& other) const { return !operator==(other); }

                    // test for exact numerical equivalence
                    constexpr bool is_exactly_the_same(const unit& other) const {
                        return base_units_ == other.base_units() && multiplier_ == other.multiplier();
                    }

                    // check if the units have the same base unit 
                    constexpr bool has_same_base(const unit& other) const { return base_units_.has_same_base(other.base_units_); }

                    // check if the units has the same base units as a unit_data object
                    constexpr bool has_same_base(unit_data base) const { return base_units_.has_same_base(base); }

                    // check if the units have the same base unit  
                    constexpr bool equivalent_non_counting(const unit& other) const {
                        return base_units_.equivalent_non_counting(other.base_units_);
                    }

                    // check if the units have the same base unit  
                    constexpr bool equivalent_non_counting(unit_data base) const {
                        return base_units_.equivalent_non_counting(base);
                    }

                    // check if the units are in some way convertible to one another
                    constexpr bool is_convertible(const unit& other) const {
                        return base_units_.equivalent_non_counting(other.base_units_);
                    } 

                    // check if the units are in some way convertible to one another
                    constexpr bool is_convertible(const unit_data& base) const {
                        return base_units_.equivalent_non_counting(base);
                    }
                    
                    // check if the unit is the default unit
                    constexpr bool is_default() const { return base_units_.empty(); }

                    // check if the multiplier is nan
                    inline bool is_nan(const unit& u) {
                        return std::isnan(u.multiplier());
                    }

                    // checks that the multiplier is finite
                    inline bool is_finite(const unit& utest) {
                        return std::isfinite(utest.multiplier());
                    }

                    // check if unit multiplier is finite
                    inline bool is_inf(const unit& utest) {
                        return std::isinf(utest.multiplier());
                    }

                    // generate a unit which is an integer power of another
                    inline constexpr unit pow(const unit& u, int power) {
                        return u.pow(power);
                    }

                    // =============================================
                    // get methods
                    // ============================================= 

                    // get the number of different base units used
                    constexpr int unit_type_count() const { return base_units_.unit_type_count(); }
                    
                    // get the base unit Multiplier
                    constexpr double multiplier() const { return multiplier_; }
                    
                    // get a rounded value of the multiplier rounded to the defined precision
                    double cround() const { return op::cround(multiplier_); }
                    
                    // get the base units
                    constexpr unit_data base_units() const { return base_units_; }


            }; // class unit     


            // verify that the units are the expected sizes
            static_assert(sizeof(unit) <= bitwidth::base_size * 2 + + sizeof(double), "Unit type is too large");


            // units declarations
            namespace defined_units {

                // some unitless numbers
                constexpr unit one;
                constexpr unit hundred = unit(100.0, one);
                constexpr unit ten = unit(10.0, one);
                constexpr unit percent(one, 0.01);
                constexpr unit infinite(unit_data(0, 0, 0, 0, 0, 0, 0), constants::infinity);
                constexpr unit neginfinite(unit_data(0, 0, 0, 0, 0, 0, 0), -constants::infinity);
                constexpr unit nan(unit_data(0, 0, 0, 0, 0, 0, 0), constants::invalid_conversion);

                // some specialized units
                constexpr unit defunit(unit_data(0, 0, 0, 0, 0, 0, 0));
                constexpr unit invalid(unit_data(nullptr), constants::invalid_conversion);
                constexpr unit error(unit_data(nullptr));
                
                
                // SI units
                namespace SI {

                    constexpr unit meter(unit_data(1, 0, 0, 0, 0, 0, 0));
                    constexpr unit m = meter;

                    constexpr unit second(unit_data(0, 1, 0, 0, 0, 0, 0));
                    constexpr unit s = second;
                    
                    constexpr unit kilogram(unit_data(0, 0, 1, 0, 0, 0, 0));
                    constexpr unit kg = kilogram;

                    constexpr unit Ampere(unit_data(0, 0, 0, 1, 0, 0, 0));
                    constexpr unit A = Ampere;

                    constexpr unit Kelvin(unit_data(0, 0, 0, 0, 1, 0, 0));
                    constexpr unit K = Kelvin;

                    constexpr unit mol(unit_data(0, 0, 0, 0, 0, 1, 0));
                    
                    constexpr unit candela(unit_data(0, 0, 0, 0, 0, 0, 1));
                    constexpr unit cd = candela;

                } // namespace SI


                // SI prefixes as units
                namespace SI_prefix {
                    
                    constexpr unit centi(one, 1e-2);
                    constexpr unit milli(one, 1e-3);
                    constexpr unit micro(one, 1e-6);
                    constexpr unit nano(one, 1e-9);
                    constexpr unit pico(one, 1e-12);
                    constexpr unit femto(one, 1e-15);
                    constexpr unit atto(one, 1e-18);
                    constexpr unit zepto(one, 1e-21);
                    constexpr unit yocto(one, 1e-24);

                    constexpr unit hecto(one, 1e2);
                    constexpr unit kilo(one, 1e3);
                    constexpr unit mega(one, 1e6);
                    constexpr unit giga(one, 1e9);
                    constexpr unit tera(one, 1e12);
                    constexpr unit peta(one, 1e15);
                    constexpr unit exa(one, 1e18);
                    constexpr unit zetta(one, 1e21);
                    constexpr unit yotta(one, 1e24);

                } // namespace SI_prefix

                    
                // derived SI units:
                namespace SI_derived { 

                    constexpr unit hertz(unit_data(0, 0, -1, 0, 0, 0, 0));
                    constexpr unit Hz = hertz;

                    constexpr unit volt(unit_data(2, 1, -3, -1, 0, 0, 0));
                    constexpr unit V = volt;

                    constexpr unit newton(unit_data(1, 1, -2, 0, 0, 0, 0));
                    constexpr unit N = newton;

                    constexpr unit Pa(unit_data(-1, 1, -2, 0, 0, 0, 0));
                    constexpr unit pascal = Pa;

                    constexpr unit joule(unit_data(2, 1, -2, 0, 0, 0, 0));
                    constexpr unit J = joule;

                    constexpr unit watt(unit_data(2, 1, -3, 0, 0, 0, 0));
                    constexpr unit W = watt;

                    constexpr unit coulomb(unit_data(0, 0, 1, 1, 0, 0, 0));
                    constexpr unit C = coulomb;

                    constexpr unit farad(unit_data(-2, -1, 4, 2, 0, 0, 0));
                    constexpr unit F = farad;

                    constexpr unit weber(unit_data(2, 1, -2, -1, 0, 0, 0));
                    constexpr unit Wb = weber;
                    
                    constexpr unit tesla(unit_data(0, 1, -2, -1, 0, 0, 0));
                    constexpr unit T = tesla;

                    constexpr unit henry(unit_data(2, 1, -2, -2, 0, 0, 0));
                    constexpr unit H = henry;

                    constexpr unit mps(SI::m / SI::s);
                    constexpr unit mpss(SI::m / SI::s.pow(2)); 

                    // distance units
                    constexpr unit km(1000.0, SI::m);
                    constexpr unit dm(0.1, SI::m);
                    constexpr unit cm(0.01, SI::m);
                    constexpr unit mm(0.001, SI::m);
                    constexpr unit um(1.e-6, SI::m);
                    constexpr unit nm(1.e-9, SI::m);

                    // volume units
                    constexpr unit L{0.001, SI::m * SI::m * SI::m};
                    constexpr unit dL{0.1, L};
                    constexpr unit cL{0.01, L};
                    constexpr unit mL{0.001, L};
                    
                    // mass units
                    constexpr unit g(0.001, SI::kg);
                    constexpr unit mg(0.001, g);

                    // time unit
                    namespace time {

                        constexpr unit ms(0.001, SI::s);
                        constexpr unit us(1.e-6, SI::s);
                        constexpr unit ns(1.e-9, SI::s);
                        constexpr unit min(60.0, SI::s);
                        constexpr unit hr(60.0, min);
                        constexpr unit day(24.0, hr);
                        constexpr unit yr(8760.0, hr);  // median calendar year;
                        constexpr unit sday{365.24 / 366.24, day};  // sidereal day
                        constexpr unit syr(365.256363004, day);  // sidereal year

                    }  // namespace time

                    constexpr unit ms = time::ms;
                    constexpr unit us = time::us; 
                    constexpr unit ns = time::ns;
                    constexpr unit min = time::min;
                    constexpr unit hr = time::hr;  
                    constexpr unit day = time::day;
                    constexpr unit yr = time::yr;
                    constexpr unit sday = time::sday; 
                    constexpr unit syr = time::syr; 


                    constexpr unit MW = SI_prefix::mega * W; 
                    constexpr unit kW = SI_prefix::kilo * W; 
                    constexpr unit mW = SI_prefix::milli * W; 
                    constexpr unit mA = SI_prefix::milli * SI::A; 
                    constexpr unit kV = SI_prefix::kilo * V; 
                    

                } // namespace SI_derived


                // // defined units to string
                // static const std::vector<unit, std::string> defined_unit_names = {
                //     // SI units
                //     {SI::m, "m"},
                //     {SI::s, "s"},
                //     {SI::kg, "kg"},
                //     {SI::A, "A"},
                //     {SI::K, "K"},
                //     {SI::mol, "mol"},
                //     {SI::cd, "cd"},
                    
                //     {SI::m * SI::m, "m^2"},
                //     {SI::m * SI::m * SI::m, "m^3"}, 
                //     {SI::s * SI::s, "s^2"},
                //     {SI::s * SI::s * SI::s, "s^3"}, 
                //     {SI::kg * SI::kg, "kg^2"},
                //     {SI::kg * SI::kg * SI::kg, "kg^3"}, 

                //     // SI derived units
                //     {SI_derived::Hz, "Hz"},
                //     {SI_derived::V, "V"},
                //     {SI_derived::N, "N"},
                //     {SI_derived::Pa, "Pa"},
                //     {SI_derived::J, "J"},
                //     {SI_derived::W, "W"},
                //     {SI_derived::C, "C"},
                //     {SI_derived::F, "F"},
                //     {SI_derived::Wb, "Wb"},
                //     {SI_derived::T, "T"},
                //     {SI_derived::H, "H"},
                //     {SI_derived::mps, "m/s"},
                //     {SI_derived::mpss, "m/s^2"},
                //     {SI_derived::MW, "MW"},
                //     {SI_derived::kW, "kW"},
                //     {SI_derived::mW, "mW"},
                //     {SI_derived::mA, "mA"},
                //     {SI_derived::kV, "kV"},


                //     {SI_derived::dm, "dm"},
                //     {SI_derived::dm * SI_derived::dm, "dm^2"},
                //     {SI_derived::dm * SI_derived::dm * SI_derived::dm, "dm^3"},
                //     {SI_derived::cm, "cm"},
                //     {SI_derived::cm * SI_derived::cm, "cm^2"},
                //     {SI_derived::cm * SI_derived::cm * SI_derived::cm, "cm^3"},
                //     {SI_derived::mm, "mm"},
                //     {SI_derived::mm * SI_derived::mm, "mm^2"},
                //     {SI_derived::mm * SI_derived::mm * SI_derived::mm, "mm^3"},
                //     {SI_derived::km, "km"},
                //     {SI_derived::km * SI_derived::km, "km^2"},
                //     {SI_derived::km * SI_derived::km * SI_derived::km, "km^3"},
                //     {SI_derived::um, "um"},
                //     {SI_derived::nm, "nm"},

                //     {SI_derived::ns, "ns"},
                //     {SI_derived::us, "us"},
                //     {SI_derived::ms, "ms"},
                //     {SI_derived::min, "min"},
                //     {SI_derived::hr, "hr"},
                //     {SI_derived::day, "day"},
                //     {SI_derived::yr, "yr"},
                    
                //     {SI_derived::g, "g"}, 
                //     {SI_derived::mg, "mg"}, 
                            
                //     {SI_derived::L, "L"},
                //     {SI_derived::dL, "dL"},
                //     {SI_derived::cL, "cL"},
                //     {SI_derived::mL, "mL"},
                    
                // };


            } // namespace defined_units


            // convertion template functions 
            namespace convert {

                // Generate a conversion factor between two units in a constexpr function, the
                // units will only convert if they have the same base unit
                template<typename UX, typename UX2>
                constexpr double quick_convert(UX start, UX2 result) { return quick_convert(1.0, start, result); }

                // Generate a conversion factor between two units in a constexpr function, the
                // units will only convert if they have the same base unit
                template<typename UX, typename UX2>
                constexpr double quick_convert(double val, const UX& start, const UX2& result) {
                    static_assert(            std::is_same<UX, unit>::value || std::is_same<UX, unit>::value,
                        "convert argument types must be unit or unit");
                    static_assert(            std::is_same<UX2, unit>::value || std::is_same<UX2, unit>::value,
                        "convert argument types must be unit or unit");
                    return (start.base_units() == result.base_units()) ?
                        val * start.multiplier() / result.multiplier() :
                            constants::invalid_conversion;
                }

                // Convert a value from one unit base to another
                template<typename UX, typename UX2>
                double convert(double val, const UX& start, const UX2& result) {
                    static_assert(std::is_same<UX, unit>::value || std::is_same<UX, unit>::value,
                        "convert argument types must be unit or unit");
                    static_assert(std::is_same<UX2, unit>::value || std::is_same<UX2, unit>::value, 
                        "convert argument types must be unit or unit");
                        
                    if (start == result) { return val; }
                    if (start.base_units() == result.base_units()) { return val * start.multiplier() / result.multiplier(); }

                    auto base_start = start.base_units();
                    auto base_result = result.base_units();

                    if (base_start.has_same_base(base_result)) { return val * start.multiplier() / result.multiplier(); }
                    if (base_start.has_same_base(base_result.inv())) { return 1.0 / (val * start.multiplier() * result.multiplier()); }
                    return constants::invalid_conversion;

                }     
                                
                // Generate a conversion factor between two units
                template<typename UX, typename UX2>
                double convert(const UX& start, const UX2& result) { return convert(1.0, start, result); }

                // static const std::unordered_map<double, char> SI_prefix_value_to_char{
                //     {0.001, 'm'},        {1.0F / 1000.0, 'm'},
                //     {1000.0, 'k'},       {1.0F / 0.001, 'k'},
                //     {1e-6, 'u'},         {1.0F / 1e6, 'u'},
                //     {0.01, 'c'},         {1.0F / 100.0, 'c'},
                //     {1000000.0, 'M'},    {1.0F / 0.000001, 'M'},
                //     {1000000000.0, 'G'}, {1.0F / 0.000000001, 'G'},
                //     {1e-9, 'n'},         {1.0F / 1e9, 'n'},
                //     {1e-12, 'p'},        {1.0F / 1e12, 'p'},
                //     {1e-15, 'f'},        {1.0F / 1e15, 'f'},
                //     {1e-18, 'a'},        {1.0F / 1e18, 'a'},
                //     {1e-21, 'z'},        {1.0F / 1e21, 'z'},
                //     {1e-24, 'y'},        {1.0F / 1e24, 'y'},
                //     {1e12, 'T'},         {1.0F / 1e-12, 'T'},
                //     {1e15, 'P'},         {1.0F / 1e-15, 'P'},
                //     {1e18, 'E'},         {1.0F / 1e-18, 'E'},
                //     {1e21, 'Z'},         {1.0F / 1e-21, 'Z'},
                //     {1e24, 'Y'},         {1.0F / 1e-24, 'Y'}
                // };

                // static const std::unordered_map<unit, char> SI_prefix_to_char {
                //     {defined_units::SI_prefix::kilo, 'k'},  
                //     {defined_units::SI_prefix::centi, 'c'}, 
                //     {defined_units::SI_prefix::milli, 'm'}, 
                //     {defined_units::SI_prefix::micro, 'u'},     
                //     {defined_units::SI_prefix::mega, 'M'},        
                //     {defined_units::SI_prefix::giga, 'G'},        
                //     {defined_units::SI_prefix::nano, 'n'},      
                //     {defined_units::SI_prefix::pico, 'p'},       
                //     {defined_units::SI_prefix::femto, 'f'},    
                //     {defined_units::SI_prefix::atto, 'a'},        
                //     {defined_units::SI_prefix::zepto, 'z'},       
                //     {defined_units::SI_prefix::yocto, 'y'},     
                //     {defined_units::SI_prefix::tera, 'T'},       
                //     {defined_units::SI_prefix::peta, 'P'},        
                //     {defined_units::SI_prefix::exa, 'E'},         
                //     {defined_units::SI_prefix::zetta, 'Z'},      
                //     {defined_units::SI_prefix::yotta, 'Y'},      
                // };

                // static const std::unordered_map<char, unit> SI_prefix_char_to_SI_prefix{
                //     {'m', defined_units::SI_prefix::milli},
                //     {'k', defined_units::SI_prefix::kilo},  
                //     {'u', defined_units::SI_prefix::micro}, 
                //     {'c', defined_units::SI_prefix::centi},      
                //     {'M', defined_units::SI_prefix::mega},        
                //     {'G', defined_units::SI_prefix::giga},        
                //     {'n', defined_units::SI_prefix::nano},      
                //     {'p', defined_units::SI_prefix::pico},       
                //     {'f', defined_units::SI_prefix::femto},    
                //     {'a', defined_units::SI_prefix::atto},        
                //     {'z', defined_units::SI_prefix::zepto},       
                //     {'y', defined_units::SI_prefix::yocto},     
                //     {'T', defined_units::SI_prefix::tera},       
                //     {'P', defined_units::SI_prefix::peta},        
                //     {'E', defined_units::SI_prefix::exa},         
                //     {'Z', defined_units::SI_prefix::zetta},      
                //     {'Y', defined_units::SI_prefix::yotta}         
                // };
                    

            } // namespace convert


        } // namespace units


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
                    constexpr measurement() noexcept {};

                    // constructor from a value and a unit 
                    constexpr measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    // constructor from a measurement
                    constexpr measurement(const measurement& other) :
                        value_(other.value()), units_(other.units()) {}

                    // constructor from a measurement
                    constexpr measurement(measurement&& val) noexcept :
                        value_(val.value()), units_(val.units()) {}
                        
                    // destructor
                    ~measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  
          
                    constexpr measurement operator*(const measurement& other) const {
                        return { value_ * other.value_, units_ * other.units_};
                    }
                    
                    constexpr measurement operator*(const units::unit& other) const {
                        return { value_, units_ * other};
                    }
                    
                    constexpr measurement operator*(double val) const {
                        return { value_ * val, units_};
                    }
                    
                    constexpr measurement operator/(const measurement& other) const {
                        return { value_ / other.value_, units_ / other.units_};
                    }
                    
                    constexpr measurement operator/(const units::unit& other) const {
                        return { value_, units_ / other};
                    }

                    constexpr measurement operator/(double val) const {
                        return { value_ / val, units_};
                    }

                    measurement operator+(const measurement& other) const {
                        return { value_ + other.value_as(units_), units_};
                    }

                    measurement operator-(const measurement& other) const {
                        return { value_ - other.value_as(units_), units_};
                    }

                    constexpr friend measurement pow(const measurement& meas, int power) {
                        return measurement{ op::power_const(meas.value_, power), meas.units_.pow(power) };
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

                    friend constexpr inline measurement operator*(double val, const units::unit& unit_base) {
                        return { val, unit_base};
                    }

                    friend constexpr inline measurement operator*(const units::unit& unit_base, double val) {
                        return { val, unit_base};
                    }

                    friend constexpr inline measurement operator/(double val, const units::unit& unit_base) {
                        return { val, unit_base.inv() };
                    }
                    
                    friend constexpr inline measurement operator/(const units::unit& unit_base, double val) {
                        return { 1.0 / val, unit_base};
                    }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return { units::convert::convert(value_, units_, newUnits), newUnits};
                    }

                    // convert a unit into its base units
                    constexpr measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::unit(units_.base_units()) };
                    }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ?
                            value_ : units::convert::convert(value_, units_, desired_units);
                    }

                    // double multiplier
                    friend constexpr inline measurement operator*(double val, const measurement& meas) {
                        return meas * val;
                    }

                    // divide measurement into a double
                    friend constexpr inline measurement operator/(double val, const measurement& meas) {
                        return { val / meas.value_, meas.units_.inv() };
                    }

                    // get the numerical component of the measurement
                    constexpr double value() const { return value_; }

                    // get the unit component of a measurement
                    constexpr units::unit units() const { return units_; }

                    // convert the measurement to a single unit
                    constexpr units::unit as_unit() const { return {value_, units_}; }

                    // print the measurement
                    void print() const {
                        std::cout << value() << " udm";  
                    }

                private:

                    // does a numerical equality check on the value accounting for tolerances
                    bool value_equality_check(double otherval) const {
                        return (value_ == otherval) ?
                            true : op::compare_round_equals(value_, otherval);
                    } 


            }; // class measurement


            // verify that the precise measurement are the expected sizes
            static_assert(sizeof(measurement) <= 2 * sizeof(double) + 2 * units::bitwidth::base_size, "precise measurement is too large");


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
                    constexpr fixed_measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    // constructor from a measurement
                    explicit constexpr fixed_measurement(const measurement& val) :
                        value_(val.value()), units_(val.units()) {}

                    // constructor from a fixed measurement
                    constexpr fixed_measurement(const fixed_measurement& val) :
                        value_(val.value()), units_(val.units()) {}

                    // constructor from a fixed measurement
                    constexpr fixed_measurement(fixed_measurement&& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    // destructor
                    ~fixed_measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

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

                    constexpr measurement operator*(const measurement& other) const {
                        return {value_ * other.value(), units_ * other.units()};
                    }

                    constexpr measurement operator*(const units::unit& other) const {
                        return {value_, units_ * other};
                    }
                    
                    constexpr fixed_measurement operator*(double val) const {
                        return {value_ * val, units_};
                    }

                    constexpr measurement operator/(const measurement& other) const {
                        return {value_ / other.value(), units_ / other.units()};
                    }

                    constexpr measurement operator/(const units::unit& other) const {
                        return {value_, units_ / other};
                    }

                    constexpr fixed_measurement operator/(double val) const {
                        return {value_ / val, units_};
                    }

                    fixed_measurement operator+(const measurement& other) const {
                        return {value_ + other.value_as(units_), units_};
                    }

                    fixed_measurement operator-(const measurement& other) const {
                        return {value_ - other.value_as(units_), units_};
                    }

                    constexpr fixed_measurement operator+(double val) const {
                        return {value_ + val, units_};
                    }

                    constexpr fixed_measurement operator-(double val) const {
                        return {value_ - val, units_};
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

                    constexpr friend fixed_measurement pow(const fixed_measurement& meas, int power) {
                        return fixed_measurement{op::power_const(meas.value_, power), meas.units_.pow(power)};
                    }

                    bool operator==(double val) const {
                        return (value_ == val) ? true : op::compare_round_equals(value_, val);
                    }

                    bool operator!=(double val) const { return !operator==(val); }

                    constexpr bool operator>(double val) const { return value_ > val; }

                    constexpr bool operator<(double val) const { return value_ < val; }

                    bool operator>=(double val) const {
                        return (value_ >= val) ? true : operator==(val);
                    }

                    bool operator<=(double val) const {
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

                    friend bool operator==(double val, const fixed_measurement& v2) {
                        return v2 == val;
                    };

                    friend bool operator!=(double val, const fixed_measurement& v2) {
                        return v2 != val;
                    };

                    friend constexpr bool operator>(double val, const fixed_measurement& v2) {
                        return val > v2.value();
                    };

                    friend constexpr bool operator<(double val, const fixed_measurement& v2) {
                        return val < v2.value();
                    };

                    friend bool operator>=(double val, const fixed_measurement& v2) {
                        return (val >= v2.value()) ? true : (v2 == val);
                    };

                    friend bool operator<=(double val, const fixed_measurement& v2) {
                        return (val <= v2.value()) ? true : (v2 == val);
                    };

                    friend constexpr inline fixed_measurement operator+(double v1, const fixed_measurement& v2) {
                        return {v1 + v2.value(), v2.units()};
                    }

                    friend constexpr inline fixed_measurement operator-(double v1, const fixed_measurement& v2) {
                        return {v1 - v2.value(), v2.units()};
                    }

                    friend constexpr inline fixed_measurement operator*(double v1, const fixed_measurement& v2) {
                        return {v1 * v2.value(), v2.units()};
                    }

                    friend constexpr inline fixed_measurement operator/(double v1, const fixed_measurement& v2) {
                        return {v1 / v2.value(), v2.units().inv()};
                    }
                    

                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // direct conversion operator
                    operator measurement() { return {value_, units_}; }

                    // set the numerical value
                    void value(double val) { value_ = val; }

                    // get the numerical value
                    constexpr double value() const { return value_; }

                    // get the unit
                    constexpr units::unit units() const { return units_; }

                    // convert the measurement to a unit
                    constexpr units::unit as_unit() const { return {value_, units_}; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? value_ : units::convert::convert(value_, units_, desired_units);
                    }

                    // convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return { units::convert::convert(value_, units_, newUnits), newUnits };
                    }

                    // convert a unit into its base units
                    constexpr measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::unit(units_.base_units()) };
                    }

                    // print the measurement
                    void print() const {
                        std::cout << value() << " udm"; 
                    }

            }; // class fixed_measurement


            // verify that the fixed measurement are the expected sizes
            static_assert(sizeof(fixed_measurement) <= 2 * sizeof(double) + 2 * units::bitwidth::base_size, "fixed measurement is too large");


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
                    constexpr uncertain_measurement() = default;
                    
                    // constructor from a value, uncertainty, and unit
                    constexpr uncertain_measurement(double val, double uncertainty_val, const units::unit& base) noexcept :
                        value_(val), uncertainty_(uncertainty_val), units_(base) {}

                    // constructpr from a value and an unit, assuming the uncertainty is 0
                    explicit constexpr uncertain_measurement(double val, const units::unit& base) noexcept :
                        value_(val), units_(base) {}

                    // constructor from a regular measurement and uncertainty value
                    explicit constexpr uncertain_measurement(const measurement& val, double uncertainty_val) noexcept : 
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
                    constexpr uncertain_measurement simple_product(const uncertain_measurement& other) const {
                        double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                        double nval = value_ * other.value_;
                        return { nval, nval * ntol, units_ * other.units() };
                    }

                    // multiply with another measurement equivalent to uncertain_measurement multiplication with 0 uncertainty
                    constexpr uncertain_measurement operator*(const measurement& other) const {
                        return { 
                            value() * other.value(),
                            other.value() * uncertainty(),
                            units_ * other.units() };
                    }

                    constexpr uncertain_measurement operator*(const units::unit& other) const {
                        return { value_, uncertainty_, units_ * other};
                    }

                    constexpr uncertain_measurement operator*(double val) const {
                        return { value_ * val, uncertainty_ * val, units_};
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
                    constexpr uncertain_measurement simple_divide(const uncertain_measurement& other) const {
                        double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                        double nval = value_ / other.value_;
                        return { nval, nval * ntol, units_ / other.units() };
                    }

                    constexpr uncertain_measurement operator/(const measurement& other) const {
                        return { 
                            static_cast<double>(value() / other.value()),
                            static_cast<double>(uncertainty() / other.value()), units_ / other.units() };
                    }

                    constexpr uncertain_measurement operator/(const units::unit& other) const {
                        return { value_, uncertainty_, units_ / other};
                    }

                    constexpr uncertain_measurement operator/(double val) const {
                        return { value_ / val, uncertainty_ / val, units_};
                    }

                    // compute a unit addition and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator+(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(units::convert::convert(other.units_, units_));
                        double ntol = std::sqrt(
                            uncertainty_ * uncertainty_ +
                            cval * cval * other.uncertainty_ * other.uncertainty_);
                        return { value_ + cval * other.value_, ntol, units_};
                    }

                    uncertain_measurement simple_add(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(units::convert::convert(other.units_, units_));
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return { value_ + cval * other.value_, ntol, units_};
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator-(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(units::convert::convert(other.units_, units_));
                        double ntol = std::sqrt(
                            uncertainty_ * uncertainty_ +
                            cval * cval * other.uncertainty_ * other.uncertainty_);
                        return { value_ - cval * other.value_, ntol, units_};
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the simple uncertainty summation method
                    uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                        auto cval = units::convert::convert(other.units_, units_);
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return { value_ - cval * other.value_, ntol, units_};
                    }

                    uncertain_measurement operator+(const measurement& other) const {
                        auto cval = other.value_as(units_);
                        return { value_ + cval, uncertainty_, units_};
                    }

                    uncertain_measurement operator-(const measurement& other) const {
                        auto cval = other.value_as(units_);
                        return { value_ - cval, uncertainty_, units_};
                    }

                    // take the measurement to some power
                    friend constexpr uncertain_measurement pow(const uncertain_measurement& meas, int power) {
                        auto new_value = op::power_const(meas.value_, power);
                        auto new_tol = ((power >= 0) ? power : -power) * new_value * meas.uncertainty_ / meas.value_;
                        return uncertain_measurement(new_value, new_tol, meas.units_.pow(power));                    
                    }

                    // comparison operators 
                    bool operator==(const measurement& other) const {
                        auto val = other.value_as(units_);
                        if (uncertainty_ == 0.0F) { return (value_ == val) ? true : op::compare_round_equals(value_, val); }
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

                    friend bool operator==(const measurement& other, const uncertain_measurement& v2) {
                        return v2 == other;
                    };

                    friend bool operator!=(const measurement& other, const uncertain_measurement& v2) {
                        return v2 != other;
                    };

                    friend constexpr bool operator>(const measurement& other, const uncertain_measurement& v2) {
                        return other.value() > v2.value();
                    };

                    friend constexpr bool operator<(const measurement& other, const uncertain_measurement& v2) {
                        return other.value() < v2.value();
                    };

                    friend bool operator>=(const measurement& other, const uncertain_measurement& v2) {
                        return (other > v2) ? true : (v2 == other);
                    };

                    friend bool operator<=(const measurement& other, const uncertain_measurement& v2) {
                        return (other < v2) ? true : (v2 == other);
                    };

                    friend inline uncertain_measurement operator+(const measurement& v1, const uncertain_measurement& v2) {
                        double cval = units::convert::convert(v2.units_, v1.units());
                        double ntol = v2.uncertainty() * cval;
                        return uncertain_measurement(
                            v1.value() + cval * v2.value(), ntol, v1.units());
                    }

                    friend inline uncertain_measurement operator-(const measurement& v1, const uncertain_measurement& v2) {
                        double cval = units::convert::convert(v2.units_, v1.units());
                        double ntol = v2.uncertainty() * cval;
                        return uncertain_measurement(v1.value() - cval * v2.value(), ntol, v1.units());
                    }

                    friend constexpr inline uncertain_measurement operator*(const measurement& v1, const uncertain_measurement& v2) {
                        return v2.operator*(v1);
                    }

                    friend constexpr inline uncertain_measurement operator/(const measurement& v1, const uncertain_measurement& v2) {
                        double ntol = v2.uncertainty() / v2.value();
                        double nval = v1.value() / v2.value();
                        return uncertain_measurement(nval, nval * ntol, v1.units() / v2.units_);
                    }

                    friend constexpr inline uncertain_measurement operator*(double v1, const uncertain_measurement& v2) {
                        return v2.operator*(v1);
                    }

                    friend constexpr inline uncertain_measurement operator*(float v1, const uncertain_measurement& v2) {
                        return v2.operator*(v1);
                    }

                    friend constexpr inline uncertain_measurement operator/(double v1, const uncertain_measurement& v2) {
                        double ntol = v2.uncertainty() / v2.value();
                        double nval = v1 / v2.value();
                        return uncertain_measurement(nval, nval * ntol, v2.units_.inv());
                    }

                    friend constexpr inline uncertain_measurement operator/(float v1, const uncertain_measurement& v2) {
                        float ntol = v2.uncertainty_ / v2.value_;
                        float nval = v1 / v2.value_;
                        return { nval, nval * ntol, v2.units_.inv() };
                    }

                    friend constexpr inline uncertain_measurement operator/(int v1, const uncertain_measurement& v2) {
                        float ntol = v2.uncertainty_ / v2.value_;
                        float nval = static_cast<float>(v1) / v2.value_;
                        return { nval, nval * ntol, v2.units_.inv() };
                    }


                    // =============================================                                                                                         
                    // get and set methods
                    // =============================================  

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
                    constexpr measurement uncertainty_measurement() const { return { uncertainty(), units_ }; }

                    // // cast operator to a measurement
                    constexpr operator measurement() const { return { value(), units_ }; }

                    // get the numerical value 
                    constexpr double value() const { return value_; }

                    // get the numerical value of the uncertainty
                    constexpr double uncertainty() const { return uncertainty_; }

                    // get the underlying units value
                    constexpr units::unit units() const { return units_; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? value_ : units::convert::convert(value_, units_, desired_units);
                    }

                    // get the numerical value of the uncertainty as a particular unit
                    double uncertainty_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ? uncertainty_ : units::convert::convert(uncertainty_, units_, desired_units);
                    }

                    // convert a unit to have a new base
                    uncertain_measurement convert_to(const units::unit& newUnits) const {
                        auto cval = units::convert::convert(units_, newUnits);
                        return { cval * value_, uncertainty_ * cval, newUnits };
                    }

                    // print the uncertain measurement
                    void print() const {
                        std::cout << value() << " ± " << uncertainty() << " udm"; 
                    }

                    // print the uncertain measurement
                    void print_as(const units::unit& desired_units) const {
                        std::cout << value_as(desired_units) << " ± " << uncertainty_as(desired_units) << " udm"; 
                    }
    
            }; // class uncertain_measurement


            // verify that the uncertain measurement are the expected sizes
            static_assert(sizeof(uncertain_measurement) <= 3 * sizeof(double) + 2 * units::bitwidth::base_size, "uncertain measurement is too large");

        } // namespace measurements
    

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

                    coordinate() = default; 


                    // =============================================
                    // set, get and print methods
                    // =============================================
                    
                    void set_coordinate(const double& x) { m_coordinate.value(x); }

                    double get_coordinate() const { return m_coordinate.value(); }

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

                    position(coordinate x, coordinate y, coordinate z) {
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

                    void set_coordinate_x(const double& x) { m_position[0].set_coordinate(x); }

                    void set_coordinate_y(const double& y) { m_position[1].set_coordinate(y); }

                    void set_coordinate_z(const double& z) { m_position[2].set_coordinate(z); }

                    std::vector<coordinate> get_coordinates() const { return m_position; }

                    coordinate get_coordinate_x() const { return m_position[0]; }

                    coordinate get_coordinate_y() const { return m_position[1]; }

                    coordinate get_coordinate_z() const { return m_position[2]; }

                    double get_magnitude() const {
                        return sqrt(pow(m_position[0].get_coordinate(), 2) +                 
                                    pow(m_position[1].get_coordinate(), 2) + 
                                    pow(m_position[2].get_coordinate(), 2));
                    }        

                    double get_distance(const std::vector<coordinate>& pos) const {        
                        return sqrt(pow(pos[0].get_coordinate() - m_position[0].get_coordinate(), 2) + 
                                    pow(pos[1].get_coordinate() - m_position[1].get_coordinate(), 2) + 
                                    pow(pos[2].get_coordinate() - m_position[2].get_coordinate(), 2)); 
                    }                    
                    
                    double get_distance(const position& pos) const {        
                        return sqrt(pow(pos.get_coordinate_x().get_coordinate() - m_position[0].get_coordinate(), 2) + 
                                    pow(pos.get_coordinate_y().get_coordinate() - m_position[1].get_coordinate(), 2) + 
                                    pow(pos.get_coordinate_z().get_coordinate() - m_position[2].get_coordinate(), 2)); 
                    }
                    
                    double get_rho() { return sqrt(pow(m_position[0].get_coordinate(), 2) + pow(m_position[1].get_coordinate(), 2)); }

                    double get_phi() const { return atan2(m_position[1].get_coordinate(), m_position[0].get_coordinate()); }     

                    double get_phi(const std::vector<coordinate>& pos) const { 
                        return atan2(pos[1].get_coordinate() - m_position[1].get_coordinate(), pos[0].get_coordinate() - m_position[0].get_coordinate()); 
                    }

                    double get_phi(const position& pos) const { 
                        return atan2(pos.get_coordinate_y().get_coordinate() - m_position[1].get_coordinate(), pos.get_coordinate_x().get_coordinate() - m_position[0].get_coordinate()); 
                    }
 
                    double get_theta() const { return acos(m_position[2].get_coordinate() / get_magnitude()); }
            
                    double get_theta(const std::vector<coordinate>& pos) const { 
                        return acos((pos[2].get_coordinate() - m_position[2].get_coordinate()) / get_distance(pos)); 
                    }

                    double get_theta(const position& pos) const { 
                        return acos((pos.get_coordinate_z().get_coordinate() - m_position[2].get_coordinate()) / get_distance(pos)); 
                    }
    
                    // std::vector<coordinate> get_direction() const {
                    //     return {cos(get_phi()), sin(get_phi()), m_position[2].get_coordinate() / get_magnitude()}};
                    // } 

                    // std::vector<coordinate> get_direction(const std::vector<coordinate>& pos1) const {
                    //     return {cos(get_phi(pos1)), sin(get_phi(pos1)), (pos1[2].get_coordinate() - m_position[2].get_coordinate()) / get_distance(pos1)};
                    // } 

                    // std::vector<coordinate> get_direction(const position& pos1) const {
                    //     return {cos(get_phi(pos1)), sin(get_phi(pos1)), (pos1.get_coordinate_z().get_coordinate() - m_position[2].get_coordinate()) / get_distance(pos1)};
                    // } 

                    // =============================================
                    // print methods
                    // =============================================

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

                    velocity(coordinate x, coordinate y, coordinate z) {
                        m_velocity.push_back(x); 
                        m_velocity.push_back(y);
                        m_velocity.push_back(z); 
                    }

                    velocity(const std::vector<coordinate>& vel) : m_velocity{vel} {}

                    velocity(const velocity& vel) : m_velocity{vel.get_velocity()} {}


                    // =============================================
                    // set, get and print methods
                    // =============================================

                    void set_velocity(const std::vector<coordinate>& vel) { m_velocity = vel; }

                    void set_velocity_x(const double& x) { m_velocity[0].set_coordinate(x); }

                    void set_velocity_y(const double& y) { m_velocity[1].set_coordinate(y);  }

                    void set_velocity_z(const double& z) { m_velocity[2].set_coordinate(z); }

                    std::vector<coordinate> get_velocity() const { return m_velocity; }

                    coordinate get_velocity_x() const { return m_velocity[0]; }

                    coordinate get_velocity_y() const { return m_velocity[1]; }

                    coordinate get_velocity_z() const { return m_velocity[2]; }
                    

                    double get_magnitude() const {
                        return sqrt(pow(m_velocity[0].get_coordinate(), 2) +                 
                                    pow(m_velocity[1].get_coordinate(), 2) + 
                                    pow(m_velocity[2].get_coordinate(), 2));
                    }          

                    double get_phi() const { return atan2(m_velocity[1].get_coordinate(), m_velocity[0].get_coordinate()); }     
                    
                    double get_theta() const { return acos(m_velocity[2].get_coordinate() / get_magnitude()); }

                    // std::vector<coordinate> get_direction() const {
                    //     return {cos(get_phi()), sin(get_phi()), cos(get_theta())};
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


        namespace elements {

            class mass {

                protected: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::fixed_measurement m_mass; 

                    position::position m_pos; 
                    
                    // bool m_gravitational_field = false;

                    // std::vector<measurements::fixed_measurement> m_gravitational_attraction = zeros(3); 


                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    mass(const double& mass, const units::unit& mass_unit, const std::vector<position::coordinate>& pos) : 
                        m_mass{mass, mass_unit}, m_pos{pos} {}

                    mass(const double& mass, const units::unit& mass_unit, const position::position& pos) : 
                        m_mass{mass, mass_unit}, m_pos{pos} {}

                    ~mass() {}


                    // // =============================================
                    // // set and get methods
                    // // =============================================

                    // void set_mass(const double& mass) { m_mass = mass; }

                    // void set_mass_um_prefix(const char* um_prefix) { m_mass_um_prefix = um_prefix; }
                    
                    // double get_mass() const { return m_mass; }

                    // const char* get_mass_um() const { return m_mass_um; }

                    // const char* get_mass_um_prefix() const { return m_mass_um_prefix; }

                
                    // // =============================================
                    // // gravitational methods
                    // // =============================================

                    // void activate_gravitational_field() { m_gravitational_field = true; }

                    // void deactivate_gravitational_field() { m_gravitational_field = false; }

                    // void reset_gravitational_attraction() { m_gravitational_attraction.clear(); }

                    // void add_gravitational_attraction(const std::vector<double>& attraction) { m_gravitational_attraction += attraction; }

                    // std::vector<double> get_gravitational_attraction() const { return m_gravitational_attraction; }

                    // std::vector<double> gravitational_attraction(const std::vector<double>& coord1) {
                    //     if (m_gravitational_field == false) {
                    //         std::cout << "Before evaluating the gravitational attraction given by this mass in these coordinates, you must activate the gravitational field." << std::endl; 
                    //         exit(-11);
                    //     }
                    //     if (coord1 == get_coordinates()) return zeros(3);
                    //     std::vector<double> appo = get_direction(coord1) * (- G * m_mass / pow(get_distance(coord1), 2));
                    //     return appo; 
                    // }


                    // // =============================================
                    // // print methods
                    // // =============================================
                    
                    // void print_mass() const { 
                    //     std::cout << "- mass = " << get_mass() << " " << get_mass_um_prefix() << get_mass_um() << std::endl;
                    // }

                    // void print_gravitational_attraction() const { 
                    //     std::cout << "- gravitational attraction: " << std::endl; 
                    //     for (auto i : get_gravitational_attraction()) std::cout << "[" << i << "]\t";
                    //     std::cout << std::endl; 
                    // }
                
            }; // class mass

        } // namespace elements

    } // namespace tools
  
} // namespace physics





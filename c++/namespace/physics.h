
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools for computational physics. 
// last updated:    23/07/2022


#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional> 
#include <limits>
#include <vector>


/* 

######## NAMESPACE PHYSICS MAP ########

namespace physics
    |
    '----> namespace tools 
    |       |
    |       '----> namespace constants     
    |       |
    |       |    
    |       '----> namespace op 
    |       |
    |       |    
    |       '----> namespace units 
    |       |       |
    |       |       '---> class unit_data
    |       |       '---> class unit 
    |       |       '---> class unit 
    |       |
    |       '----> namespace vector_algebra
    |       |       |
    |       |       '--->
    |       |       '--->    
    |       |       '--->
    |       |       '--->
    |       |
    |       '----> namespace measurements
    |               |
    |               '--->
    |               '--->    
    |               '--->
    |               '--->
    |
*/



namespace physics {

    namespace tools {

        namespace constants {

            constexpr int32_t max_neg(uint32_t n_bits) { return -(int32_t(1U << (n_bits - 1))); }
            constexpr double infinity = std::numeric_limits<double>::infinity();
            constexpr double invalid_conversion = std::numeric_limits<double>::signaling_NaN();

            constexpr double pi = 3.14159265358979323846;

        } // namespace constants

        static_assert(
            std::numeric_limits<double>::has_signaling_NaN, 
                "nan is used to signify invalid values");
        
        static_assert(
            std::numeric_limits<double>::has_infinity, 
                "nan is used to signify invalid values");

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
                    
                    // construct from powers
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
                    inline bool isnan(const unit& u) {
                        return std::isnan(u.multiplier());
                    }

                    // checks that the multiplier is finite
                    inline bool isfinite(const unit& utest) {
                        return std::isfinite(utest.multiplier());
                    }

                    // check if unit multiplier is finite
                    inline bool isinf(const unit& utest) {
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
                    static_assert(
                        std::is_same<UX, unit>::value || std::is_same<UX, unit>::value,
                        "convert argument types must be unit or unit");
                    static_assert(
                        std::is_same<UX2, unit>::value || std::is_same<UX2, unit>::value,
                        "convert argument types must be unit or unit");
                    return (start.base_units() == result.base_units()) ?
                        val * start.multiplier() / result.multiplier() :
                            constants::invalid_conversion;
                }

                // Generate a conversion factor between two units
                template<typename UX, typename UX2>
                double convert(const UX& start, const UX2& result) { return convert(1.0, start, result); }

                // Convert a value from one unit base to another
                template<typename UX, typename UX2>
                double convert(double val, const UX& start, const UX2& result) {
                    static_assert( 
                        std::is_same<UX, unit>::value || std::is_same<UX, unit>::value,
                        "convert argument types must be unit or unit");
                    static_assert( 
                        std::is_same<UX2, unit>::value || std::is_same<UX2, unit>::value, 
                        "convert argument types must be unit or unit");
                        
                    if (start == result) { return val; }
                    if (start.base_units() == result.base_units()) { return val * start.multiplier() / result.multiplier(); }

                    auto base_start = start.base_units();
                    auto base_result = result.base_units();

                    if (base_start.has_same_base(base_result)) { return val * start.multiplier() / result.multiplier(); }
                    if (base_start.has_same_base(base_result.inv())) { return 1.0 / (val * start.multiplier() * result.multiplier()); }
                    return constants::invalid_conversion;

                }    

            } // namespace convert


            // units declarations
            namespace defined_units {

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


                    // // define some specialized units
                    constexpr unit defunit(unit_data(0, 0, 0, 0, 0, 0, 0));
                    constexpr unit invalid(unit_data(nullptr), constants::invalid_conversion);
                    constexpr unit error(unit_data(nullptr));

                    // Define some unitless numbers
                    constexpr unit one;
                    constexpr unit hundred = unit(100.0, one);
                    constexpr unit ten = unit(10.0, one);
                    constexpr unit percent(one, 0.01);
                    constexpr unit infinite(unit_data(0, 0, 0, 0, 0, 0, 0), constants::infinity);
                    constexpr unit neginfinite(unit_data(0, 0, 0, 0, 0, 0, 0), -constants::infinity);
                    constexpr unit nan(unit_data(0, 0, 0, 0, 0, 0, 0), constants::invalid_conversion);

                // SI prefixes as units
                namespace SI_prefix {

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

                    constexpr unit kibi(one, 1024);
                    constexpr unit mebi = kibi * kibi;
                    constexpr unit gibi = mebi * kibi;
                    constexpr unit tebi = gibi * kibi;
                    constexpr unit pebi = tebi * kibi;
                    constexpr unit exbi = pebi * kibi;
                    constexpr unit zebi = exbi * kibi;
                    constexpr unit yobi = zebi * kibi;

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

                    constexpr unit mps(m / s);


                    // Distance units
                    constexpr unit cm(m, 0.01);
                    constexpr unit km(m, 1000.0);
                    constexpr unit mm(m, 0.001);
                    constexpr unit nm(m, 1e-9);

                    // Volume units
                    constexpr unit L{m * m * m, 0.001};
                    constexpr unit mL{L, 0.001};
                    
                    // mass units
                    constexpr unit g(kg, 0.001);
                    constexpr unit mg(g, 0.001);

                    namespace time {

                        // Time unit
                        constexpr unit min(s, 60.0);
                        constexpr unit ms(s, 0.001);
                        constexpr unit ns(s, 1e-9);
                        constexpr unit hr(min, 60.0);
                        constexpr unit h(min, 60.0);
                        constexpr unit day(hr, 24.0);
                        constexpr unit week(day, 7.0);
                        constexpr unit fortnight(day, 14);
                        constexpr unit yr(hr, 8760.0);  // median calendar year;
                        constexpr unit year = yr;  // standard year for SI
                        constexpr unit sday{day, 365.24 / 366.24};  // sidereal day
                        constexpr unit syr(day, 365.256363004);  // sidereal year

                    }  // namespace time

                    constexpr unit min = time::min;
                    constexpr unit ms = time::ms;
                    constexpr unit ns = time::ns;
                    constexpr unit hr = time::hr;
                    constexpr unit h = time::h;
                    constexpr unit yr = time::yr;
                    constexpr unit day = time::day;

                }

            } // namespace 

        } // namespace units


        namespace measurements {

            // class using precise units and double precision
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

                    // constructor from a value and a fixed unit 
                    constexpr measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    // constructor from implicit conversion from a lower precision measurement
                    constexpr measurement(const measurement& other) :
                        value_(other.value()), units_(other.units()) {}

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

                private:

                    // does a numerical equality check on the value accounting for tolerances
                    bool value_equality_check(double otherval) const {
                        return (value_ == otherval) ?
                            true : op::compare_round_equals(value_, otherval);
                    } 


            }; // class measurement


            // verify that the precise measurement are the expected sizes
            static_assert(
                sizeof(measurement) <= 2 * sizeof(double) + 2 * units::bitwidth::base_size,
                    "precise measurement is too large");


            // class using precise units and a value
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

                    constexpr fixed_measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    explicit constexpr fixed_measurement(const measurement& val) :
                        value_(val.value()), units_(val.units()) {}

                    constexpr fixed_measurement(const fixed_measurement& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

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
                        return (value_ == val) ?
                            true : op::compare_round_equals(value_, val);
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
                        return operator==(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator!=(const fixed_measurement& val) const {
                        return operator!=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator==(const measurement& val) const {
                        return operator==(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator!=(const measurement& val) const {
                        return operator!=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator>(const measurement& val) const {
                        return operator>(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator<(const measurement& val) const {
                        return operator<(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator>=(const measurement& val) const {
                        return operator>=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator<=(const measurement& val) const {
                        return operator<=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
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

                    constexpr double value() const { return value_; }

                    constexpr units::unit units() const { return units_; }

                    // convert the measurement to a single unit
                    constexpr units::unit as_unit() const { return {value_, units_}; }

                    // Get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ?
                            value_ :
                            units::convert::convert(value_, units_, desired_units);
                    }

                    // Convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return { units::convert::convert(value_, units_, newUnits), newUnits };
                    }

                    // Convert a unit into its base units
                    constexpr measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::unit(units_.base_units()) };
                    }


                }; // class fixed_measurement


            // verify that the fixed measurement are the expected sizes
            static_assert(
                sizeof(fixed_measurement) <= sizeof(double) + 2 * units::bitwidth::base_size,
                    "fixed measurement is too large");


            // class using precise units and an uncertain value
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

                    constexpr uncertain_measurement() = default;
                    
                    // construct from a single precision value, uncertainty, and unit
                    constexpr uncertain_measurement(double val, double uncertainty_val, const units::unit& base) noexcept :
                        value_(val), uncertainty_(uncertainty_val), units_(base) {}

                    // construct from a single precision value, and unit assume uncertainty is 0
                    explicit constexpr uncertain_measurement(double val, const units::unit& base) noexcept :
                        value_(val), units_(base) {}

                    // construct from a regular measurement and uncertainty value
                    explicit constexpr uncertain_measurement(const measurement& val, double uncertainty_val) noexcept : 
                        value_(val.value()), uncertainty_(uncertainty_val), units_(val.units()) {}

                    // construct from a regular measurement and an uncertainty measurement
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
                            static_cast<double>(uncertainty() / other.value()),
                            units_ / other.units() };
                    }

                    constexpr uncertain_measurement operator/(const units::unit& other) const {
                        return { value_, uncertainty_, units_ / other};
                    }

                    constexpr uncertain_measurement operator/(double val) const {
                        return { value_ / val, uncertainty_ / val, units_};
                    }

                    // compute a unit addition and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator+(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(convert(other.units_, units_));
                        double ntol = std::sqrt(
                            uncertainty_ * uncertainty_ +
                            cval * cval * other.uncertainty_ * other.uncertainty_);
                        return { value_ + cval * other.value_, ntol, units_};
                    }

                    uncertain_measurement simple_add(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(convert(other.units_, units_));
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return { value_ + cval * other.value_, ntol, units_};
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator-(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(convert(other.units_, units_));
                        double ntol = std::sqrt(
                            uncertainty_ * uncertainty_ +
                            cval * cval * other.uncertainty_ * other.uncertainty_);
                        return { value_ - cval * other.value_, ntol, units_};
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the simple uncertainty summation method
                    uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                        auto cval = convert::convert(other.units_, units_);
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
                        return uncertain_measurement{ new_value, new_tol, meas.units_.pow(power) };
                    }

                    // comparison operators 
                    bool operator==(const measurement& other) const {
                        auto val = other.value_as(units_);
                        if (uncertainty_ == 0.0F) { return (value_ == val) ? true :
                                                    op::compare_round_equals(value_, val);
                        }
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
                        double cval = convert::convert(v2.units_, v1.units());
                        double ntol = v2.uncertainty() * cval;
                        return uncertain_measurement(
                            v1.value() + cval * v2.value(), ntol, v1.units());
                    }

                    friend inline uncertain_measurement operator-(const measurement& v1, const uncertain_measurement& v2) {
                        double cval = convert::convert(v2.units_, v1.units());
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

                    // Convert a unit to have a new base
                    uncertain_measurement convert_to(const units::unit& newUnits) const {
                        auto cval = convert::convert(units_, newUnits);
                        return { cval * value_, uncertainty_ * cval, newUnits};
                    }

                    // get the underlying units value
                    constexpr units::unit units() const { return units_; }

                    // get the numerical value as a particular unit type
                    double value_as(units::unit desired_units) const {
                        return (units_ == desired_units) ? value_ :  units::convert::convert(value_, units_, desired_units);
                    }
                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return value_as(units::unit_cast(desired_units));
                    }
                    // get the numerical value of the uncertainty as a particular unit
                    double uncertainty_as(units::unit desired_units) const {
                        return (units_ == desired_units) ? uncertainty_ : units::convert::convert(uncertainty_, units_, desired_units);
                    }

                    double uncertainty_as(const units::unit& desired_units) const {
                        return uncertainty_as(units::unit_cast(desired_units));
                    }


                    // get the base value with no units
                    constexpr double value() const { return value_; }
                    
                    // get the uncertainty with no units
                    constexpr double uncertainty() const { return uncertainty_; }

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

                    // get the fractional uncertainty with no units
                    constexpr double fractional_uncertainty() const {
                        return uncertainty() / ((value_ >= 0.0F) ? value_ : -value_);
                    }

                    // get the uncertainty as a separate measurement
                    constexpr measurement uncertainty_measurement() const { return { uncertainty(), units_}; }

                    // cast operator to a measurement
                    constexpr operator measurement() const { return { value(), units_}; }


            }; // class uncertain_measurement


            // verify that the uncertain measurement are the expected sizes
            static_assert(
                sizeof(uncertain_measurement) <= 2 * sizeof(double) + 2 * units::bitwidth::base_size,
                    "uncertain measurement is too large");


        } // namespace measurements

    
    } // namespace tools
  

} // namespace physics

            


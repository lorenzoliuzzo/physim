
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools for computational physics. 
// last updated:    14/07/2022


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
    |       '----> namespace units 
    |       |       |
    |       |       '---> class unit_data
    |       |       '---> class base 
    |       |       '---> enum class prefix_enum
    |       |       '---> class prefix 
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

            // round the multiplier to the expected level of precision
            inline float cround(float val) {
                std::uint32_t bits;
                std::memcpy(&bits, &val, sizeof(bits));
                bits += 8UL;
                bits &= 0xFFFFFFF0UL;
                std::memcpy(&val, &bits, sizeof(bits));
                return val;
            }

            // round a value to the expected level of precision of a double
            inline double cround_precise(double val) {
                std::uint64_t bits;
                std::memcpy(&bits, &val, sizeof(bits));
                bits += 0x800ULL;
                bits &= 0xFFFFFFFFFFFFF000ULL;
                std::memcpy(&val, &bits, sizeof(bits));
                return val;
            }

            // rounding compare for equality on floats
            inline bool compare_round_equals(float val1, float val2) {
                static constexpr float half_precision{5e-7F};
                auto v1 = val1 - val2;
                if (v1 == 0.0F || std::fpclassify(v1) == FP_SUBNORMAL) return true;
                auto c1 = cround(val1);
                auto c2 = cround(val2);
                return (c1 == c2) || 
                    (cround(val2 * (1.0F + half_precision)) == c1) ||
                    (cround(val2 * (1.0F - half_precision)) == c1) ||
                    (cround(val1 * (1.0F + half_precision)) == c2) ||
                    (cround(val1 * (1.0F - half_precision)) == c2);
            }

            // rounding compare for equality on double
            inline bool compare_round_equals_precise(double val1, double val2) {
                static constexpr double half_precise_precision{5e-13};
                auto v1 = val1 - val2;
                if (v1 == 0.0 || std::fpclassify(v1) == FP_SUBNORMAL) { return true; }
                auto c1 = cround_precise(val1);
                auto c2 = cround_precise(val2);
                return (c1 == c2) ||
                    (cround_precise(val2 * (1.0 + half_precise_precision)) == c1) ||
                    (cround_precise(val2 * (1.0 - half_precise_precision)) == c1) ||
                    (cround_precise(val1 * (1.0 + half_precise_precision)) == c2) ||
                    (cround_precise(val1 * (1.0 - half_precise_precision)) == c2);
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
                            -candela_};
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
                                                                candela_ / power) : unit_data(nullptr);
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

                    // Check if the units have the same base unit 
                    constexpr bool has_same_base(const unit_data& other) const {
                        return equivalent_non_counting(other) && mole_ == other.mole_;
                    }

                    // Check equivalence for non-counting base units
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

            };    

            // class defining a basic unit module with float precision on the multiplier
            class unit {

                private:

                    // =============================================
                    // class members
                    // ============================================= 

                    friend class precise_unit;

                    unit_data base_units_{0, 0, 0, 0, 0, 0, 0};
                    
                    float multiplier_{1.0};


                public:
                    
                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // default constructor
                    constexpr unit() noexcept {};
                    
                    // construct unit from base unit 
                    explicit constexpr unit(const unit_data& base_unit) : base_units_(base_unit) {};

                    // construct unit from base unit and a multiplier
                    constexpr unit(const unit_data& base_unit, double mult) : base_units_(base_unit), multiplier_(static_cast<float>(mult)) {};

                    // construct unit from base unit and a multiplier
                    constexpr explicit unit(const unit_data& base_unit, float mult) : base_units_(base_unit), multiplier_(mult) {};  
                    
                    // construct unit from multiplier and a base unit 
                    constexpr unit(double mult, const unit& other) : unit(other.base_units_, mult * other.multiplier()) {};
                    

                    // =============================================
                    // operations
                    // ============================================= 
    
                    // unit multiplication
                    constexpr unit operator*(const unit& other) const { return { base_units_ * other.base_units_, multiplier() * other.multiplier() }; }

                    // division operator
                    constexpr unit operator/(const unit& other) const { return { base_units_ / other.base_units_, multiplier() / other.multiplier() }; }

                    // invert the unit 
                    constexpr unit inv() const { return { base_units_.inv(), 1.0 / multiplier() }; }

                    // take a unit to an integral power
                    constexpr unit pow(int power) const { return unit{base_units_.pow(power), op::power_const(multiplier_, power) }; }

                    
                    // =============================================
                    // checks
                    // ============================================= 

                    // test for unit equivalence to within nominal numerical tolerance (6 decimal digits)
                    bool operator==(const unit& other) const {
                        if (base_units_ != other.base_units_) { return false; }
                        if (multiplier_ == other.multiplier_) { return true; }
                        return op::compare_round_equals(multiplier_, other.multiplier_);
                    }

                    bool operator!=(const unit& other) const { return !operator == (other); }

                    // test for exact numerical equivalence
                    constexpr bool is_exactly_the_same(const unit& other) const { return base_units_ == other.base_units_ && multiplier_ == other.multiplier_; }

                    // check if the units have the same base unit 
                    constexpr bool has_same_base(const unit& other) const { return base_units_.has_same_base(other.base_units_); }

                    constexpr bool has_same_base(const unit_data& base) const { return base_units_.has_same_base(base); }
    
                    // check if the units have the same base unit 
                    constexpr bool equivalent_non_counting(const unit& other) const { return base_units_.equivalent_non_counting(other.base_units_); }

                    constexpr bool equivalent_non_counting(const unit_data& base) const { return base_units_.equivalent_non_counting(base); }
    
                    // check if the units are in some way convertible to one another

                    constexpr bool is_convertible(const unit& other) const { return base_units_.equivalent_non_counting(other.base_units_); }

                    constexpr bool is_convertible(const unit_data& base) const { return base_units_.equivalent_non_counting(base); }


                    // =============================================
                    // get methods
                    // =============================================     
                    
                    // get the number of different base units used
                    constexpr int unit_type_count() const { return base_units_.unit_type_count(); }    
                    
                    // extract the base unit Multiplier
                    constexpr double multiplier() const { return static_cast<double>(multiplier_); }
                    
                    // extract the base unit Multiplier as a single precision float
                    constexpr float multiplier_f() const { return multiplier_; }
                    
                    // generate a rounded version of the multiplier
                    float cround() const { return op::cround(multiplier_); }

                    // get the unit_data
                    constexpr unit_data base_units() const { return base_units_; }

            };

            // class defining a basic unit module with double precision on the multiplier
            class precise_unit {
                
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
                    constexpr precise_unit() noexcept {};

                    // construct from base_unit 
                    explicit constexpr precise_unit(const unit_data& base_unit) noexcept : 
                        base_units_(base_unit) {};
                    
                    // construct from base_unit and multiplier
                    constexpr precise_unit(const unit_data& base_unit, double mult) noexcept :
                        base_units_(base_unit), multiplier_(mult) {}
                    
                    // copy constructor from a less precise unit
                    explicit constexpr precise_unit(const unit& other) noexcept : 
                        base_units_(other.base_units_), multiplier_(other.multiplier()) {};

                    // copy constructor with a multiplier
                    constexpr precise_unit(const precise_unit& other, double mult) noexcept :
                        precise_unit(other.base_units_, mult * other.multiplier_) {};

                    // constructor from unit with a double 
                    constexpr precise_unit(const unit& other, double mult) noexcept :
                        precise_unit(other.base_units(), mult * other.multiplier()) {};

                    // Constructor from a double with an unit
                    constexpr precise_unit(double mult, const unit& other) noexcept :
                        precise_unit(other.base_units(), mult * other.multiplier()) {};


                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a unit
                    constexpr precise_unit inv() const { return { base_units_.inv(), 1.0 / multiplier() }; }

                    // multiply with another unit
                    constexpr precise_unit operator*(const precise_unit& other) const {
                        return {  base_units_ * other.base_units_, 
                            multiplier() * other.multiplier() };
                    }

                    // multiply with a lower precision unit
                    constexpr precise_unit operator*(const unit& other) const {
                        return {  base_units_ * other.base_units_,
                            multiplier() * other.multiplier() };
                    }

                    // division operator
                    constexpr precise_unit operator/(const precise_unit& other) const {
                        return {  base_units_ / other.base_units_,
                            multiplier() / other.multiplier() };
                    }

                    // divide by a less precise unit
                    constexpr precise_unit operator/(const unit& other) const {
                        return {  base_units_ / other.base_units_, 
                            multiplier() / other.multiplier() };
                    }

                    // take a unit to a power
                    constexpr precise_unit pow(int power) const {
                        return {  base_units_.pow(power), op::power_const(multiplier_, power) };
                    }


                    // =============================================
                    // checks
                    // ============================================= 

                    // equality operator
                    bool operator==(const precise_unit& other) const {
                        if (base_units_ != other.base_units_ ) { return false; }
                        if (multiplier_ == other.multiplier_) { return true; }
                        return op::compare_round_equals_precise(multiplier_, other.multiplier_);
                    }

                    bool operator!=(const precise_unit& other) const { return !operator==(other); }

                    // test for exact numerical equivalence
                    constexpr bool is_exactly_the_same(const unit& other) const {
                        return base_units_ == other.base_units() && multiplier_ == other.multiplier();
                    }

                    // Check if the units have the same base unit (i.e. they measure the same
                    // thing)
                    constexpr bool has_same_base(const precise_unit& other) const { return base_units_.has_same_base(other.base_units_); }

                    // Check if the units have the same base unit (i.e. they measure the same
                    // thing)
                    constexpr bool has_same_base(const unit& other) const { return base_units_.has_same_base(other.base_units_); }

                    // Check if the units has the same base units as a unit_data object
                    constexpr bool has_same_base(unit_data base) const { return base_units_.has_same_base(base); }

                    // Check rounded equality with another unit
                    bool operator==(const unit& other) const {
                        if (base_units_ != other.base_units_) { return false; }
                        if (multiplier_ == other.multiplier()) { return true; }
                        return op::compare_round_equals(static_cast<float>(multiplier_), other.multiplier_);
                    }

                    bool operator!=(const unit& other) const { return !operator==(other); }

                    // Check if the units have the same base unit (i.e. they measure the same
                    // thing)
                    
                    constexpr bool equivalent_non_counting(const precise_unit& other) const {
                        return base_units_.equivalent_non_counting(other.base_units_);
                    }
                    
                    constexpr bool equivalent_non_counting(const unit& other) const {
                        return base_units_.equivalent_non_counting(other.base_units_);
                    }
                    
                    constexpr bool equivalent_non_counting(unit_data base) const {
                        return base_units_.equivalent_non_counting(base);
                    }

                    // Check if the units are in some way convertible to one another
                    
                    constexpr bool is_convertible(const precise_unit& other) const {
                        return base_units_.equivalent_non_counting(other.base_units_);
                    }
                    
                    constexpr bool is_convertible(const unit& other) const {
                        return base_units_.equivalent_non_counting(other.base_units_);
                    }
                    
                    constexpr bool is_convertible(const unit_data& base) const {
                        return base_units_.equivalent_non_counting(base);
                    }
                    
                    // Check unit equality (base units equal and equivalent multipliers to
                    // specified precision
                    friend bool operator==(const unit& val1, const precise_unit& val2) { return (val2 == val1); }
                    
                    // get the number of different base units used
                    constexpr int unit_type_count() const { return base_units_.unit_type_count(); }

                    // check if the unit is the default unit
                    constexpr bool is_default() const { return base_units_.empty(); }
                    
                    // Extract the base unit Multiplier
                    constexpr double multiplier() const { return multiplier_; }
                    
                    // Extract the base unit Multiplier as a single precision float
                    constexpr float multiplier_f() const { return static_cast<float>(multiplier_); }
                    
                    // Generate a rounded value of the multiplier rounded to the defined
                    // precision
                    double cround() const { return op::cround_precise(multiplier_); }
                    
                    // get the base units
                    constexpr unit_data base_units() const { return base_units_; }

            }; // class precise_unit     
            
            // check if a unit down cast is lossless
            inline constexpr bool is_unit_cast_lossless(const precise_unit& val) {
                return val.multiplier() == static_cast<double>(static_cast<float>(val.multiplier()));
            }

            // downcast a precise unit to the less precise version
            constexpr unit unit_cast(const precise_unit& val) {
                return { val.base_units(), val.multiplier() };
            }

            // downcast a unit to the less precise version
            constexpr unit unit_cast(const unit& val) { return val; }

            // Check if the multiplier is nan
            inline bool isnan(const precise_unit& u) {
                return std::isnan(u.multiplier());
            }

            // Check if the multiplier is nan
            inline bool isnan(const unit& u) {
                return std::isnan(u.multiplier_f());
            }

            //checks that the multiplier is finite
            inline bool isfinite(const precise_unit& utest) {
                return std::isfinite(utest.multiplier());
            }

            // check if the unit multiplier is finite
            inline bool isfinite(const unit& utest) {
                return std::isfinite(utest.multiplier_f());
            }

            // check if unit multiplier is finite
            inline bool isinf(const precise_unit& utest) {
                return std::isinf(utest.multiplier());
            }

            // check if unit multiplier is infinite
            inline bool isinf(const unit& utest) {
                return std::isinf(utest.multiplier_f());
            }

            // generate a unit which is an integer power of another
            inline constexpr unit pow(const unit& u, int power) {
                return u.pow(power);
            }

            // generate a precise unit which is an integer power of another
            inline constexpr precise_unit pow(const precise_unit& u, int power) {
                return u.pow(power);
            }

            // Verify that the units are the expected sizes
            static_assert(
                sizeof(unit_data) == bitwidth::base_size, 
                    "Unit data is too large");

            static_assert(
                sizeof(unit) <= bitwidth::base_size * 2,
                    "Unit type is too large");

            static_assert(
                sizeof(precise_unit) <= bitwidth::base_size * 2 + sizeof(double),
                    "precise unit type is too large");


            namespace precise {

                // base units
                constexpr precise_unit meter(unit_data(1, 0, 0, 0, 0, 0, 0));
                constexpr precise_unit m = meter;

                constexpr precise_unit second(unit_data(0, 1, 0, 0, 0, 0, 0));
                constexpr precise_unit s = second;
                
                constexpr precise_unit kilogram(unit_data(0, 0, 1, 0, 0, 0, 0));
                constexpr precise_unit kg = kilogram;

                constexpr precise_unit Ampere(unit_data(0, 0, 0, 1, 0, 0, 0));
                constexpr precise_unit A = Ampere;

                constexpr precise_unit Kelvin(unit_data(0, 0, 0, 0, 1, 0, 0));
                constexpr precise_unit K = Kelvin;

                constexpr precise_unit mol(unit_data(0, 0, 0, 0, 0, 1, 0));
                
                constexpr precise_unit candela(unit_data(0, 0, 0, 0, 0, 0, 1));
                constexpr precise_unit cd = candela;

                // // define some specialized units
                constexpr precise_unit defunit(unit_data(0, 0, 0, 0, 0, 0, 0));
                constexpr precise_unit invalid(unit_data(nullptr), constants::invalid_conversion);
                constexpr precise_unit error(unit_data(nullptr));

                // Define some unitless numbers
                constexpr precise_unit one;
                constexpr precise_unit hundred = precise_unit(one, 100.0);
                constexpr precise_unit ten = precise_unit(one, 10.0);
                constexpr precise_unit percent(one, 0.01);

                constexpr precise_unit infinite(unit_data(0, 0, 0, 0, 0, 0, 0), constants::infinity);
                constexpr precise_unit neginfinite(unit_data(0, 0, 0, 0, 0, 0, 0), -constants::infinity);
                constexpr precise_unit nan(unit_data(0, 0, 0, 0, 0, 0, 0), constants::invalid_conversion);

                // SI prefixes as units
                constexpr precise_unit milli(one, 1e-3);
                constexpr precise_unit micro(one, 1e-6);
                constexpr precise_unit nano(one, 1e-9);
                constexpr precise_unit pico(one, 1e-12);
                constexpr precise_unit femto(one, 1e-15);
                constexpr precise_unit atto(one, 1e-18);
                constexpr precise_unit zepto(one, 1e-21);
                constexpr precise_unit yocto(one, 1e-24);

                constexpr precise_unit hecto(one, 1e2);
                constexpr precise_unit kilo(one, 1e3);
                constexpr precise_unit mega(one, 1e6);
                constexpr precise_unit giga(one, 1e9);
                constexpr precise_unit tera(one, 1e12);
                constexpr precise_unit peta(one, 1e15);
                constexpr precise_unit exa(one, 1e18);
                constexpr precise_unit zetta(one, 1e21);
                constexpr precise_unit yotta(one, 1e24);

                constexpr precise_unit kibi(one, 1024);
                constexpr precise_unit mebi = kibi * kibi;
                constexpr precise_unit gibi = mebi * kibi;
                constexpr precise_unit tebi = gibi * kibi;
                constexpr precise_unit pebi = tebi * kibi;
                constexpr precise_unit exbi = pebi * kibi;
                constexpr precise_unit zebi = exbi * kibi;
                constexpr precise_unit yobi = zebi * kibi;

                // Derived SI units:
                constexpr precise_unit hertz(unit_data(0, 0, -1, 0, 0, 0, 0));
                constexpr precise_unit Hz = hertz;

                constexpr precise_unit volt(unit_data(2, 1, -3, -1, 0, 0, 0));
                constexpr precise_unit V = volt;

                constexpr precise_unit newton(unit_data(1, 1, -2, 0, 0, 0, 0));
                constexpr precise_unit N = newton;

                constexpr precise_unit Pa(unit_data(-1, 1, -2, 0, 0, 0, 0));
                constexpr precise_unit pascal = Pa;

                constexpr precise_unit joule(unit_data(2, 1, -2, 0, 0, 0, 0));
                constexpr precise_unit J = joule;

                constexpr precise_unit watt(unit_data(2, 1, -3, 0, 0, 0, 0));
                constexpr precise_unit W = watt;

                constexpr precise_unit coulomb(unit_data(0, 0, 1, 1, 0, 0, 0));
                constexpr precise_unit C = coulomb;

                constexpr precise_unit farad(unit_data(-2, -1, 4, 2, 0, 0, 0));
                constexpr precise_unit F = farad;

                constexpr precise_unit weber(unit_data(2, 1, -2, -1, 0, 0, 0));
                constexpr precise_unit Wb = weber;
                
                constexpr precise_unit tesla(unit_data(0, 1, -2, -1, 0, 0, 0));
                constexpr precise_unit T = tesla;

                constexpr precise_unit henry(unit_data(2, 1, -2, -2, 0, 0, 0));
                constexpr precise_unit H = henry;

                constexpr precise_unit mps(m / s);


                // Distance units
                constexpr precise_unit cm(m, 0.01);
                constexpr precise_unit km(m, 1000.0);
                constexpr precise_unit mm(m, 0.001);
                constexpr precise_unit nm(m, 1e-9);

                // Volume units
                constexpr precise_unit L{m * m * m, 0.001};
                constexpr precise_unit mL{L, 0.001};
                
                // mass units
                constexpr precise_unit g(kg, 0.001);
                constexpr precise_unit mg(g, 0.001);

                namespace time {

                    // Time unit
                    constexpr precise_unit min(s, 60.0);
                    constexpr precise_unit ms(s, 0.001);
                    constexpr precise_unit ns(s, 1e-9);
                    constexpr precise_unit hr(min, 60.0);
                    constexpr precise_unit h(min, 60.0);
                    constexpr precise_unit day(hr, 24.0);
                    constexpr precise_unit week(day, 7.0);
                    constexpr precise_unit fortnight(day, 14);
                    constexpr precise_unit yr(hr, 8760.0);  // median calendar year;
                    constexpr precise_unit year = yr;  // standard year for SI
                    constexpr precise_unit sday{day, 365.24 / 366.24};  // sidereal day
                    constexpr precise_unit syr(day, 365.256363004);  // sidereal year

                }  // namespace time

                constexpr precise_unit min = time::min;
                constexpr precise_unit ms = time::ms;
                constexpr precise_unit ns = time::ns;
                constexpr precise_unit hr = time::hr;
                constexpr precise_unit h = time::h;
                constexpr precise_unit yr = time::yr;
                constexpr precise_unit day = time::day;
 
            } // namespace precise

            // Generate a conversion factor between two units in a constexpr function, the
            // units will only convert if they have the same base unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(UX start, UX2 result) { return quick_convert(1.0, start, result); }

            // Generate a conversion factor between two units in a constexpr function, the
            // units will only convert if they have the same base unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(double val, const UX& start, const UX2& result) {
                static_assert(
                    std::is_same<UX, unit>::value || std::is_same<UX, precise_unit>::value,
                    "convert argument types must be unit or precise_unit");
                static_assert(
                    std::is_same<UX2, unit>::value ||
                        std::is_same<UX2, precise_unit>::value,
                    "convert argument types must be unit or precise_unit");
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
                    std::is_same<UX, unit>::value || std::is_same<UX, precise_unit>::value,
                    "convert argument types must be unit or precise_unit");
                static_assert( 
                    std::is_same<UX2, unit>::value || std::is_same<UX2, precise_unit>::value, 
                    "convert argument types must be unit or precise_unit");
                    
                if (start == result) { return val; }
                if (start.base_units() == result.base_units()) { return val * start.multiplier() / result.multiplier(); }

                auto base_start = start.base_units();
                auto base_result = result.base_units();

                if (base_start.has_same_base(base_result)) { return val * start.multiplier() / result.multiplier(); }
                if (base_start.has_same_base(base_result.inv())) { return 1.0 / (val * start.multiplier() * result.multiplier()); }
                return constants::invalid_conversion;

            }

            // template<typename UX, typename UX2>
            // inline double otherUsefulConversions(double val, const UX& start, const UX2& result) {
            //     if (start.has_same_base(N) && result.has_same_base(kg)) {
            //         // weight to mass
            //         return val * start.multiplier() / standard_gravity / result.multiplier();
            //     }
            //     if (start.has_same_base(kg) && result.has_same_base(N)) {
            //         // mass to weight
            //         return val * start.multiplier() * standard_gravity / result.multiplier();
            //     }
            //     return constants::invalid_conversion;
            // }

        } // namespace units

            
        namespace measurements {

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

                    // Default constructor
                    constexpr measurement() noexcept {}

                    // construct from a value and unit
                    constexpr measurement(double val, const units::unit& base) : value_(val), units_(base) {}


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

                    measurement operator%(const measurement& other) const {
                        return { fmod(value_, other.value_as(units_)), units_};
                    }

                    measurement operator%(double val) const {
                        return { fmod(value_, val), units_};
                    }

                    measurement operator+(const measurement& other) const {
                        return { value_ + other.value_as(units_), units_};
                    }

                    measurement operator-(const measurement& other) const {
                        return { value_ - other.value_as(units_), units_};
                    }

                    // double multiplier
                    friend constexpr inline measurement operator*(double val, const measurement& meas) {
                        return meas * val;
                    }

                    // divide measurement into a double
                    friend constexpr inline measurement operator/(double val, const measurement& meas) {
                        return { val / meas.value_, meas.units_.inv() };
                    }

                    friend constexpr measurement pow(const measurement& meas, int power) {
                        return {  op::power_const(meas.value_, power), meas.units_.pow(power) };
                    }
                    
                    // Equality operator
                    bool operator==(const measurement& other) const {
                        auto val = other.value_as(units_);
                        return (value_ == val) ?
                            true : op::compare_round_equals(static_cast<float>(value_), static_cast<float>(val));
                    }

                    bool operator>(const measurement& other) const {
                        return value_ > other.value_as(units_);
                    }

                    bool operator<(const measurement& other) const {
                        return value_ < other.value_as(units_);
                    }

                    bool operator>=(const measurement& other) const {
                        auto val = other.value_as(units_);
                        return (value_ >= val) ?
                            true : op::compare_round_equals(static_cast<float>(value_), static_cast<float>(val));
                    }

                    bool operator<=(const measurement& other) const {
                        auto val = other.value_as(units_);
                        return (value_ <= val) ?
                            true : op::compare_round_equals(static_cast<float>(value_), static_cast<float>(val));
                    }

                    // Not equal operator
                    bool operator!=(const measurement& other) const {
                        return !operator==(other);
                    }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // Convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return { convert(value_, units_, newUnits), newUnits};
                    }

                    // Convert a unit into its base units
                    constexpr measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::unit(units_.base_units()) };
                    }

                    friend constexpr inline measurement operator*(double val, const units::unit& unit_base) { return { val, unit_base}; }
                    
                    friend constexpr inline measurement operator*(const units::unit& unit_base, double val) { return { val, unit_base}; }

                    friend constexpr inline measurement operator/(double val, const units::unit& unit_base) { return { val, unit_base.inv() }; }
                    
                    friend constexpr inline measurement operator/(const units::unit& unit_base, double val) { return { 1.0 / val, unit_base}; }


                    // =============================================                                                                                         
                    // get methods
                    // =============================================  

                    // get the base value with no units
                    constexpr double value() const { return value_; }

                    // extract the current units from the measurement
                    constexpr units::unit units() const { return units_; }

                    // convert the measurement to a single unit
                    constexpr units::unit as_unit() const { return { value_, units_}; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_unit) const {
                        return (units_ == desired_unit) ?
                            value_ : convert(value_, units_, desired_unit);
                    }

                    // get the numerical value as a particular unit type
                    double value_as(units::precise_unit desired_units) const {
                        return value_as(unit_cast(desired_units));
                    }

            }; // class measurement

            // Verify that the measurement are the expected sizes
            static_assert(
                sizeof(measurement) <= 2 * units::bitwidth::base_size + sizeof(double),
                    "Measurement class is too large");

            class fixed_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0};  

                    const units::unit units_;  // a fixed unit of measurement


                public:

                    // =============================================                                                                                         
                    // constructors and destructor
                    // =============================================  

                    // construct from a value and unit
                    constexpr fixed_measurement(double val, const units::unit& base) :
                        value_(val), units_(base) {}

                    // construct from a regular measurement
                    explicit constexpr fixed_measurement(const measurement& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    // copy constructor
                    constexpr fixed_measurement(const fixed_measurement& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    // move constructor
                    constexpr fixed_measurement(fixed_measurement&& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    // destructor
                    ~fixed_measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    // assignment operator
                    fixed_measurement& operator=(const measurement& val) {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }
                    // assignment operator treat it the same as a measurement
                    fixed_measurement& operator=(const fixed_measurement& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    // assignment operator treat it the same as a measurement
                    fixed_measurement& operator=(fixed_measurement&& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    constexpr measurement operator*(const measurement& other) const {
                        return { value_ * other.value(), units_ * other.units() };
                    }
                    constexpr measurement operator*(const units::unit& other) const {
                        return { value_, units_ * other};
                    }
                    constexpr fixed_measurement operator*(double val) const {
                        return { value_ * val, units_};
                    }
                    constexpr measurement operator/(const measurement& other) const {
                        return { value_ / other.value(), units_ / other.units() };
                    }
                    constexpr measurement operator/(const units::unit& other) const {
                        return { value_, units_ / other};
                    }
                    constexpr fixed_measurement operator/(double val) const {
                        return { value_ / val, units_};
                    }

                    fixed_measurement operator+(const measurement& other) const {
                        return { value_ + other.value_as(units_), units_};
                    }
                    fixed_measurement operator-(const measurement& other) const {
                        return { value_ - other.value_as(units_), units_};
                    }

                    constexpr fixed_measurement operator+(double val) const {
                        return { value_ + val, units_};
                    }
                    constexpr fixed_measurement operator-(double val) const {
                        return { value_ - val, units_};
                    }

                    // take the measurement to some power
                    constexpr friend fixed_measurement pow(const fixed_measurement& meas, int power) {
                        return {  op::power_const(meas.value_, power), meas.units_.pow(power) };
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

                    // comparison operators
                    bool operator==(double val) const {
                        return (value_ == val) ?
                            true : op::compare_round_equals(static_cast<float>(value_), static_cast<float>(val));
                    }

                    bool operator!=(double val) const { return !operator==(val); }
                    
                    constexpr bool operator>(double val) const { return value_ > val; }
                    
                    constexpr bool operator<(double val) const { return value_ < val; }
                    
                    bool operator>=(double val) const { return (value_ > val) ? true : operator==(val); }
                    
                    bool operator<=(double val) const { return value_ < val ? true : operator==(val); }

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

                    friend bool operator==(double val, const fixed_measurement& v2) { return v2 == val; }
                    
                    friend bool operator!=(double val, const fixed_measurement& v2) { return v2 != val; }
                    
                    friend constexpr bool operator>(double val, const fixed_measurement& v2) { return val > v2.value(); }

                    friend constexpr bool operator<(double val, const fixed_measurement& v2) { return val < v2.value(); }

                    friend bool operator>=(double val, const fixed_measurement& v2) { return (val > v2.value()) ? true : (v2 == val); }

                    friend bool operator<=(double val, const fixed_measurement& v2) { return (val < v2.value()) ? true : (v2 == val); }

                    // friend operators for math operators
                    friend constexpr inline fixed_measurement operator+(double v1, const fixed_measurement& v2) {
                        return { v1 + v2.value(), v2.units() };
                    }
                    friend constexpr inline fixed_measurement operator-(double v1, const fixed_measurement& v2) {
                        return { v1 - v2.value(), v2.units() };
                    }
                    friend constexpr inline fixed_measurement operator*(double v1, const fixed_measurement& v2) {
                        return { v1 * v2.value(), v2.units() };
                    }
                    friend constexpr inline fixed_measurement operator/(double v1, const fixed_measurement& v2) {
                        return { v1 / v2.value(), v2.units().inv() };
                    }


                    // =============================================                                                                                         
                    // get methods
                    // =============================================  

                    // direct conversion operator
                    operator measurement() { return { value_, units_}; }

                    // Convert a unit to have a new base
                    fixed_measurement convert_to(const units::unit& newUnits) const {
                        return { units::convert(value_, units_, newUnits), newUnits};
                    }

                    // get the underlying units value
                    constexpr units::unit units() const { return units_; }

                    // convert the measurement to a single unit
                    constexpr units::unit as_unit() const { return { value_, units_}; }

                    // get the base value with no units
                    constexpr double value() const { return value_; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (units_ == desired_units) ?
                            value_ :
                            units::convert(static_cast<double>(value_), units_, desired_units);
                    }

                    // get the numerical value as a particular unit type
                    double value_as(units::precise_unit desired_units) const {
                        return value_as(units::unit_cast(desired_units));
                    }

            }; // class fixed_measurement       

            // Verify that the fixed measurement are the expected sizes
            static_assert(
                sizeof(fixed_measurement) <= sizeof(double) + 2 * units::bitwidth::base_size,
                    "fixed measurement is too large");


            class uncertain_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0}; 

                    double uncertainty_{0.0};

                    units::unit units_;  // fixed unit of measurement        
                

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
                    uncertain_measurement
                        simple_subtract(const uncertain_measurement& other) const {
                        auto cval = convert(other.units_, units_);
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

                    // Not equal operator
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

                    // Not equal operator
                    bool operator!=(const uncertain_measurement& other) const {
                        return !operator==(other);
                    }

                    // operator for measurement =
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

                    // friend operators for math operators
                    friend inline uncertain_measurement operator+(const measurement& v1, const uncertain_measurement& v2) {
                        double cval = convert(v2.units_, v1.units());
                        double ntol = v2.uncertainty() * cval;
                        return uncertain_measurement(
                            v1.value() + cval * v2.value(), ntol, v1.units());
                    }

                    friend inline uncertain_measurement operator-(const measurement& v1, const uncertain_measurement& v2) {
                        double cval = convert(v2.units_, v1.units());
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
                        auto cval = convert(units_, newUnits);
                        return { cval * value_, uncertainty_ * cval, newUnits};
                    }

                    // Get the underlying units value
                    constexpr units::unit units() const { return units_; }

                    // Get the numerical value as a particular unit type
                    double value_as(units::unit desired_units) const {
                        return (units_ == desired_units) ? value_ :  units::convert(value_, units_, desired_units);
                    }
                    // Get the numerical value as a particular unit type
                    double value_as(const units::precise_unit& desired_units) const {
                        return value_as(units::unit_cast(desired_units));
                    }
                    // Get the numerical value of the uncertainty as a particular unit
                    double uncertainty_as(units::unit desired_units) const {
                        return (units_ == desired_units) ? uncertainty_ : units::convert(uncertainty_, units_, desired_units);
                    }

                    double uncertainty_as(const units::precise_unit& desired_units) const {
                        return uncertainty_as(units::unit_cast(desired_units));
                    }


                    // get the base value with no units
                    constexpr double value() const { return value_; }
                    
                    // get the uncertainty with no units
                    constexpr double uncertainty() const { return uncertainty_; }

                    // get the base value with no units as a single precision double
                    constexpr double value_f() const { return value_; }

                    // get the uncertainty with no units as a single precision double
                    constexpr double uncertainty_f() const { return uncertainty_; }

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

            // Verify that the uncertain measurement are the expected sizes
            static_assert(
                sizeof(uncertain_measurement) <= 2 * sizeof(double) + 2 * units::bitwidth::base_size,
                    "uncertain measurement is too large");


            // Class using precise units and double precision
            class precise_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0};

                    units::precise_unit units_;


                public:

                    // =============================================                                                                                         
                    // constructors
                    // =============================================  

                    // Default constructor
                    constexpr precise_measurement() noexcept {};

                    constexpr precise_measurement(double val, const units::precise_unit& base) :
                        value_(val), units_(base) {}

                    // implicit conversion from a lower precision measurement
                    constexpr precise_measurement(const measurement& other) :
                        value_(other.value()), units_(other.units()) {}


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    constexpr precise_measurement operator*(const precise_measurement& other) const {
                        return { value_ * other.value_, units_ * other.units_};
                    }
                    constexpr precise_measurement operator*(const units::precise_unit& other) const {
                        return { value_, units_ * other};
                    }
                    constexpr precise_measurement operator*(double val) const {
                        return { value_ * val, units_};
                    }
                    constexpr precise_measurement operator/(const precise_measurement& other) const {
                        return { value_ / other.value_, units_ / other.units_};
                    }
                    constexpr precise_measurement operator/(const units::precise_unit& other) const {
                        return { value_, units_ / other};
                    }
                    constexpr precise_measurement operator/(double val) const {
                        return { value_ / val, units_};
                    }

                    precise_measurement operator+(const precise_measurement& other) const {
                        return { value_ + other.value_as(units_), units_};
                    }
                    precise_measurement operator-(const precise_measurement& other) const {
                        return { value_ - other.value_as(units_), units_};
                    }

                    // take the measurement to some power
                    constexpr friend precise_measurement pow(const precise_measurement& meas, int power) {
                        return precise_measurement{ op::power_const(meas.value_, power), meas.units_.pow(power) };
                    }

                    // Equality operator
                    bool operator==(const precise_measurement& other) const {
                        return valueEqualityCheck((units_ == other.units()) ? other.value() : other.value_as(units_));
                    }
                    // Not equal operator
                    bool operator!=(const precise_measurement& other) const {
                        return !valueEqualityCheck((units_ == other.units()) ? other.value() : other.value_as(units_));
                    }

                    bool operator>(const precise_measurement& other) const {
                        return value_ > other.value_as(units_);
                    }
                    bool operator<(const precise_measurement& other) const {
                        return value_ < other.value_as(units_);
                    }
                    bool operator>=(const precise_measurement& other) const {
                        double val = other.value_as(units_);
                        return (value_ > val) ? true : valueEqualityCheck(val);
                    }
                    bool operator<=(const precise_measurement& other) const {
                        double val = other.value_as(units_);
                        return (value_ < val) ? true : valueEqualityCheck(val);
                    }

                    friend constexpr inline precise_measurement operator*(double val, const units::precise_unit& unit_base) {
                        return { val, unit_base};
                    }

                    friend constexpr inline precise_measurement operator*(const units::precise_unit& unit_base, double val) {
                        return { val, unit_base};
                    }

                    friend constexpr inline precise_measurement operator/(double val, const units::precise_unit& unit_base) {
                        return { val, unit_base.inv() };
                    }
                    
                    friend constexpr inline precise_measurement operator/(const units::precise_unit& unit_base, double val) {
                        return { 1.0 / val, unit_base};
                    }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // Convert a unit to have a new base
                    precise_measurement convert_to(const units::precise_unit& newUnits) const {
                        return { units::convert(value_, units_, newUnits), newUnits};
                    }

                    // Convert a unit into its base units
                    constexpr precise_measurement convert_to_base() const {
                        return { value_ * units_.multiplier(), units::precise_unit(units_.base_units()) };
                    }

                    // Get the numerical value as a particular unit type
                    double value_as(const units::precise_unit& desired_units) const {
                        return (units_ == desired_units) ?
                            value_ : units::convert(value_, units_, desired_units);
                    }
                    // double multiplier
                    friend constexpr inline precise_measurement operator*(double val, const precise_measurement& meas) {
                        return meas * val;
                    }
                    // divide measurement into a double
                    friend constexpr inline precise_measurement operator/(double val, const precise_measurement& meas) {
                        return { val / meas.value_, meas.units_.inv() };
                    }

                    // getter for the numerical component of the measurement
                    constexpr double value() const { return value_; }

                    // getter for the unit component of a measurement
                    constexpr units::precise_unit units() const { return units_; }

                    // convert the measurement to a single unit
                    constexpr units::precise_unit as_unit() const { return {value_, units::unit_cast(units_)}; }

                private:

                    // does a numerical equality check on the value accounting for tolerances
                    bool valueEqualityCheck(double otherval) const {
                        return (value_ == otherval) ?
                            true : op::compare_round_equals_precise(value_, otherval);
                    } 

            }; // class precise_measurement

            // Verify that the precise measurement are the expected sizes
            static_assert(
                sizeof(precise_measurement) <= 2 * sizeof(double) + 2 * units::bitwidth::base_size,
                    "precise measurement is too large");


            /// Class using precise units and double precision
            class fixed_precise_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{0.0};  

                    const units::precise_unit units_;  


                public:

                    // =============================================                                                                                         
                    // constructors
                    // =============================================  

                    constexpr fixed_precise_measurement(double val, const units::precise_unit& base) :
                        value_(val), units_(base) {}

                    explicit constexpr fixed_precise_measurement(const precise_measurement& val) :
                        value_(val.value()), units_(val.units()) {}

                    constexpr fixed_precise_measurement(const fixed_precise_measurement& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    constexpr fixed_precise_measurement(fixed_precise_measurement&& val) noexcept :
                        value_(val.value()), units_(val.units()) {}

                    /// destructor
                    ~fixed_precise_measurement() = default;

                    /// assign from a precise_measurement do the conversion and assign the value
                    fixed_precise_measurement& operator=(const precise_measurement& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    /// assign from another fixed_precise_measurement treat like a
                    /// precise_measurement
                    fixed_precise_measurement& operator=(const fixed_precise_measurement& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    fixed_precise_measurement& operator=(fixed_precise_measurement&& val) noexcept {
                        value_ = (units_ == val.units()) ? val.value() : val.value_as(units_);
                        return *this;
                    }

                    /// Assignment from double,  allow direct numerical assignment since the
                    /// units are fixes and known at construction time
                    fixed_precise_measurement& operator=(double val) noexcept {
                        value_ = val;
                        return *this;
                    }
                    /// direct conversion operator
                    // NOLINTNEXTLINE(google-explicit-constructor)
                    operator precise_measurement() { return {value_, units_}; }

                    constexpr double value() const { return value_; }
                    constexpr precise_unit units() const { return units_; }

                    // convert the measurement to a single unit
                    constexpr precise_unit as_unit() const { return {value_, units_}; }

                    /// Get the numerical value as a particular unit type
                    double value_as(const precise_unit& desired_units) const {
                        return (units_ == desired_units) ?
                            value_ :
                            units::convert(value_, units_, desired_units);
                    }

                    constexpr precise_measurement operator*(const precise_measurement& other) const {
                        return {value_ * other.value(), units_ * other.units()};
                    }
                    constexpr precise_measurement operator*(const precise_unit& other) const {
                        return {value_, units_ * other};
                    }
                    
                    constexpr fixed_precise_measurement operator*(double val) const {
                        return {value_ * val, units_};
                    }
                    constexpr precise_measurement operator/(const precise_measurement& other) const {
                        return {value_ / other.value(), units_ / other.units()};
                    }
                    constexpr precise_measurement operator/(const precise_unit& other) const {
                        return {value_, units_ / other};
                    }
                    constexpr fixed_precise_measurement operator/(double val) const {
                        return {value_ / val, units_};
                    }

                    fixed_precise_measurement operator+(const precise_measurement& other) const {
                        return {value_ + other.value_as(units_), units_};
                    }
                    fixed_precise_measurement operator-(const precise_measurement& other) const {
                        return {value_ - other.value_as(units_), units_};
                    }
                    /// Add a double assuming the same units
                    constexpr fixed_precise_measurement operator+(double val) const {
                        return {value_ + val, units_};
                    }
                    /// Subtract a double assuming the same units
                    constexpr fixed_precise_measurement operator-(double val) const {
                        return {value_ - val, units_};
                    }

                    fixed_precise_measurement& operator+=(double val) {
                        value_ += val;
                        return *this;
                    }
                    fixed_precise_measurement& operator-=(double val) {
                        value_ -= val;
                        return *this;
                    }
                    fixed_precise_measurement& operator*=(double val) {
                        value_ *= val;
                        return *this;
                    }
                    fixed_precise_measurement& operator/=(double val) {
                        value_ /= val;
                        return *this;
                    }

                    /// take the measurement to some power
                    constexpr friend fixed_precise_measurement
                        pow(const fixed_precise_measurement& meas, int power) {
                        return fixed_precise_measurement{
                            detail::power_const(meas.value_, power), meas.units_.pow(power)};
                    }

                    /// Convert a unit to have a new base
                    precise_measurement convert_to(const precise_unit& newUnits) const {
                        return {units::convert(value_, units_, newUnits), newUnits};
                    }

                    /// Convert a unit into its base units
                    constexpr precise_measurement convert_to_base() const {
                        return {
                            value_ * units_.multiplier(), precise_unit(units_.base_units())};
                    }

                    /// comparison operators
                    bool operator==(double val) const {
                        return (value_ == val) ?
                            true :
                            detail::compare_round_equals_precise(value_, val);
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

                    /// Equality operator
                    bool operator==(const fixed_precise_measurement& val) const {
                        return operator==(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }
                    /// Not equal operator
                    bool operator!=(const fixed_precise_measurement& val) const {
                        return operator!=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    /// Equality operator
                    bool operator==(const precise_measurement& val) const {
                        return operator==(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }
                    /// Not equal operator
                    bool operator!=(const precise_measurement& val) const {
                        return operator!=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    bool operator>(const precise_measurement& val) const {
                        return operator>(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }
                    bool operator<(const precise_measurement& val) const {
                        return operator<(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }
                    bool operator>=(const precise_measurement& val) const {
                        return operator>=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }
                    bool operator<=(const precise_measurement& val) const {
                        return operator<=(
                            (units_ == val.units()) ? val.value() : val.value_as(units_));
                    }

                    friend bool operator==(double val, const fixed_precise_measurement& v2) {
                        return v2 == val;
                    };
                    friend bool operator!=(double val, const fixed_precise_measurement& v2) {
                        return v2 != val;
                    };
                    friend constexpr bool operator>(double val, const fixed_precise_measurement& v2) {
                        return val > v2.value();
                    };
                    friend constexpr bool operator<(double val, const fixed_precise_measurement& v2) {
                        return val < v2.value();
                    };
                    friend bool operator>=(double val, const fixed_precise_measurement& v2) {
                        return (val >= v2.value()) ? true : (v2 == val);
                    };
                    friend bool operator<=(double val, const fixed_precise_measurement& v2) {
                        return (val <= v2.value()) ? true : (v2 == val);
                    };

                    /// friend operators for math operators
                    friend constexpr inline fixed_precise_measurement operator+(double v1, const fixed_precise_measurement& v2) {
                        return {v1 + v2.value(), v2.units()};
                    }
                    friend constexpr inline fixed_precise_measurement operator-(double v1, const fixed_precise_measurement& v2) {
                        return {v1 - v2.value(), v2.units()};
                    }
                    friend constexpr inline fixed_precise_measurement operator*(double v1, const fixed_precise_measurement& v2) {
                        return {v1 * v2.value(), v2.units()};
                    }
                    friend constexpr inline fixed_precise_measurement operator/(double v1, const fixed_precise_measurement& v2) {
                        return {v1 / v2.value(), v2.units().inv()};
                    }


                }; // namespace fixed_prefix_measurement

        } // namespace measurements 







        // namespace op {

        //     template<typename X>
        //     X numericalRoot(X value, int power) {
        //         switch (power) {
        //             case 0: return X{1.0};
        //             case 1: return value;
        //             case -1: return X{1.0} / value;
        //             case 2: 
        //                 if (value < X{0.0}) { return constants::invalid_conversion; }
        //                 return std::sqrt(value);
        //             case -2:
        //                 if (value < X{0.0}) { return constants::invalid_conversion; }
        //                 return std::sqrt(X{1.0} / value);
        //             case 3:
        //                 return std::cbrt(value);
        //             case -3:
        //                 return std::cbrt(X{1.0} / value);
        //             case 4:
        //                 if (value < X{0.0}) { return constants::invalid_conversion; }
        //                 return std::sqrt(std::sqrt(value));
        //             case -4:
        //                 if (value < X{0.0}) { return constants::invalid_conversion; }
        //                 return std::sqrt(std::sqrt(X{1.0} / value));
        //             default:
        //                 if (value < X{0.0} && power % 2 == 0) { return constants::invalid_conversion; }
        //                 return std::pow(value, X{1.0} / static_cast<X>(power));
        //         }
        //     }

        //     units::unit root(const units::unit& un, int power) {
        //         if (power == 0) { return units::precise::one; }
        //         if (un.multiplier() < 0.0 && power % 2 == 0) { return error; }
        //         return unit{ un.base_units().root(power), numericalRoot(un.multiplier(), power) };
        //     }

        //     precise_unit root(const precise_unit& un, int power) {
        //         if (power == 0) { return precise::one; }
        //         if (un.multiplier() < 0.0 && power % 2 == 0) { return precise::invalid; }
        //         return precise_unit{ un.base_units().root(power), numericalRoot(un.multiplier(), power) };
        //     }

        //     measurement root(const measurement& meas, int power) {
        //         return { numericalRoot(meas.value(), power), root(meas.units(), power) };
        //     }

        //     fixed_measurement root(const fixed_measurement& fm, int power) {
        //         return { numericalRoot(fm.value(), power), root(fm.units(), power) };
        //     }

        //     uncertain_measurement root(const uncertain_measurement& um, int power) {
        //         auto new_value = numericalRoot(um.value(), power);
        //         auto new_tol = new_value * um.uncertainty() /
        //             (static_cast<double>((power >= 0) ? power : -power) * um.value());
        //         return uncertain_measurement(new_value, new_tol, root(um.units(), power));
        //     }

        //     precise_measurement root(const precise_measurement& pm, int power) {
        //         return { numericalRoot(pm.value(), power), root(pm.units(), power) };
        //     }

        //     fixed_precise_measurement root(const fixed_precise_measurement& fpm, int power) {
        //         return { numericalRoot(fpm.value(), power), root(fpm.units(), power) };
        //     }


        // } // namespace op











// constexpr std::array<std::pair<unit, const char*>, 100> defined_unit_names{
//         {{m, "m"},
//          {m * m, "m^2"},
//          {m * m * m, "m^3"},
//          {kg, "kg"},
//          {mol, "mol"},
//          {unit(0.1, m), "dm"},  // don't want to use deci in most conversion
//                                 // contexts but for meter and liter it is fine
//          {unit(0.1, L), "dL"},
//          {A, "A"},
//          {s, "s"},
//          {cd, "cd"},
//          {K, "K"},
//          {N, "N"},
//          {Pa, "Pa"},
//          {J, "J"},
//          {C, "C"},
//          {F, "F"},
//          {Wb, "Wb"},
//          {T, "T"},
//          {H, "H"},
//          {bar, "bar"},
//          {deg, "deg"},
//          {rad, "rad"},
//          {min, "min"},
//          {ms, "ms"},
//          {ns, "ns"},
//          {h, "h"},
//          {unit_cast(precise::time::day), "day"},
//          {unit_cast(precise::time::week), "week"},
//          {unit_cast(precise::time::yr), "yr"},
//          {unit_cast(precise::time::syr), "syr"},
//          {cm, "cm"},
//          {km, "km"},
//          {km * km, "km^2"},
//          {mm, "mm"},
//          {nm, "nm"},
//          {percent, "%"},
//          {error, "ERROR"},
//          {rpm, "rpm"},
//          {W, "W"},
//          {MW, "MW"},
//          {kW, "kW"},
//          {mW, "mW"},
//          {puMW, "puMW"},
//          {mA, "mA"},
//          {kV, "kV"},
//          {unit_cast(precise::pressure::atm), "atm"},
//          {L, "L"},
//          {unit_cast(precise::mL), "mL"},
//          {g, "g"},
//          {mg, "mg"},
//         {"1", precise::one},
//         {"one", precise::one},
//         {"inf", precise::infinite},
//         {"-inf", precise::neginfinite},
//         {"infinite", precise::infinite},
//         {"nan", precise::nan},
//         {"NaN", precise::nan},
//         {"deci", precise_unit(0.1, precise::one)},
//         {"0.1", precise_unit(0.1, precise::one)},
//         {".1", precise_unit(0.1, precise::one)},
//         {"deci", precise_unit(0.1, precise::one)},
//         {"0.01", precise_unit(0.01, precise::one)},
//         {".01", precise_unit(0.01, precise::one)},
//         {"centi", precise_unit(0.01, precise::one)},
//         {"0.001", precise::milli},
//         {".001", precise::milli},
//         {"milli", precise::milli},
//         {"1e-3", precise::milli},
//         {"1e-6", precise::micro},
//         {"micro", precise::micro},
//         {"1e-9", precise::nano},
//         {"1e-12", precise::pico},
//         {"1e-15", precise::femto},
//         {"1e-18", precise::atto},
//         {"nano", precise::nano},
//         {"pico", precise::pico},
//         {"femto", precise::femto},
//         {"atto", precise::atto},
//         {"10", precise_unit(10.0, precise::one)},
//         {"100", precise_unit(100, precise::one)},
//         {"hundred", precise_unit(100, precise::one)},
//         {"fifty", precise_unit(50, precise::one)},
//         {"1000", precise::kilo},
//         {"thousand", precise::kilo},
//         {"1000000", precise::mega},
//         {"million", precise::mega},
//         {"1000000000", precise::giga},
//         {"billion", precise::giga},
//         {"trillion", precise::tera},
//         {"quadrillion", precise::peta},
//         {"1e3", precise::kilo},
//         {"1e6", precise::mega},
//         {"1e9", precise::giga},
//         {"1e12", precise::tera},
//         {"1e15", precise::peta},
//         {"1e18", precise::exa},
//         {"kilo", precise::kilo},
//         {"mega", precise::mega},
//         {"giga", precise::giga},
//         {"tera", precise::tera},
//         {"peta", precise::peta},
//         {"exa", precise::exa},
//         {"%", precise::percent},
//         {"percent", precise::percent},
//         {"percentage", precise::percent},
//         {"permille", precise::milli},
//         {"m", precise::m},
//         {"meter", precise::m},
//         {"squaremeter", precise::m.pow(2)},  // to simplify some things later in the chain
//         {"cubicmeter", precise::m.pow(3)},
//         {"micron", precise::micro* precise::m},
//         {"cubiccentimeter", precise::cm.pow(3)},
//         {"m/s^2", precise::m / precise::s.pow(2)},
//         {"kg/m^3", precise::kg / precise::m.pow(3)},
//         {"kg", precise::kg},
//         {"KG", precise::kg},
//         {"kilogram", precise::kg},
//         {"mol", precise::mol},
//         {"grampercent", precise_unit(10.0, precise::g / precise::L)},
//         {"V", precise::V},
//         {"volt", precise::V},
//         {"W", precise::W},
//         {"W/m^2", precise::W / precise::m.pow(2)},
//         {"watt", precise::W},
//         {"kW", precise::electrical::kW},
//         {"kilowatt", precise::electrical::kW},
//         {"MW", precise::MW},
//         {"megawatt", precise::MW},
//         {"second", precise::s},
//         {"cd", precise::cd},
//         {"candela", precise::cd},
//         {"K", precise::K},
//         {"kelvin", precise::K},
//         {"degKelvin", precise::K},
//         {"degK", precise::K},
//         {"N", precise::N},
//         {"Ns", precise::N * precise::s},  // this would not pass through to the
//                                          // separation functions
//         {"Nm", precise::N * precise::m},  // this would not pass through to the
//         {"newton", precise::N},
//         {"Pa", precise::Pa},
//         {"pa", precise::Pa},
//         {"pascal", precise::Pa},
//         {"J", precise::J},
//         {"joule", precise::J},
//         {"Joule", precise::J},
//         {"C", precise::C},
//         {"coulomb", precise::C},
//         {"faraday", precise::other::faraday},
//         {"ohm", precise::ohm},
//         {"Ohm", precise::ohm},
//         {"Wb", precise::Wb},
//         {"weber", precise::Wb},
//         {"T", precise::T},
//         {"tesla", precise::T},
//         {"H", precise::H},
//         {"henry", precise::H},
//         {"min", precise::min},
//         {"ms", precise::ms},
//         {"millisecond", precise::ms},
//         {"hr", precise::hr},
//         {"HR", precise::hr},
//         {"h", precise::hr},
//         {"hour", precise::hr},
//         {"day", precise::time::day},
//         {"week", precise::time::week},
//         {"wk", precise::time::week},
//         {"y", precise::time::year},
//         {"yr", precise::time::yr},  // this one gets 365 days exactly
//         {"year", precise::time::year},  // year SI Definition 365 days
//         {"decade", precise::ten* precise::time::aj},  // year
//         {"century", precise::hundred* precise::time::aj},  // year
//         {"millennia", precise::kilo* precise::time::ag},  // year
//         {"millennium", precise::kilo* precise::time::ag},  // year
//         {"rad", precise::rad},
//     }
// }













        namespace vector_algebra {

            using std::vector;
            
            // =============================================                                                                                         
            // Utilities
            // =============================================  

            vector<double> zeros(const int& n) { return vector<double>(n, 0.); }

            vector<vector<double>> zeros(const int& n_rows, const int& n_cols) { return vector<vector<double>>(n_rows, zeros(n_cols)); }

            void print(const vector<vector<double>>& mat) {
                for (int i{}; i < mat.size(); i++) {
                    for (int j{}; j < mat.front().size(); j++) {
                        std::cout << "[" << mat[i][j] << "]\t";
                    }
                    std::cout << std::endl; 
                }   
            }

            // =============================================                                                                                         
            // Sum
            // =============================================  

            // sum of a vector and a scalar
            vector<double> operator+(const vector<double>& vec1, const double& value) {
                vector<double> vec = zeros(vec1.size());
                for (unsigned int i{}; i < vec1.size(); i++) vec[i] = vec1[i] + value;
                return vec;
            }

            // sum of a scalar and a vector  
            vector<double> operator+(const double& value, const vector<double>& vec1) {
                vector<double> vec = zeros(vec1.size());
                for (unsigned int i{}; i < vec1.size(); i++) vec[i] = value + vec1[i];
                return vec;
            }

            // sum of two vectors
            vector<double> operator+(const vector<double>& vec1, const vector<double>& vec2) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] + vec2[i];
                return vec;
            }

            // sum of a multidimentional vector and a scalar
            vector<vector<double>> operator+(const vector<vector<double>>& mat1, const double& value) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = mat1[i][j] + value;
                    }
                }
                return mat;
            }

            // sum of a scalar and a multidimentional vector 
            vector<vector<double>> operator+(const double& value, const vector<vector<double>>& mat1) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = value + mat1[i][j];
                    }
                }
                return mat;
            }

            // sum of two multidimentional vectors
            vector<vector<double>> operator+(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = mat1[i][j] + mat2[i][j];
                    }
                }
                return mat;
            }

            // increase vector with a scalar
            void operator+=(vector<double> vec1, const double& value) {
                for (int i{}; i < vec1.size(); i++) vec1[i] += value;
            }

            // increase vector with a vector
            void operator+=(vector<double> vec1, const vector<double>& vec2) {
                for (int i{}; i < vec1.size(); i++) vec1[i] += vec2[i];
            }

            // increase multidimentional vector with a scalar
            void operator+=(vector<vector<double>> mat1, const double& value) {
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1[i].size(); j++) {
                        mat1[i][j] += value;
                    }
                }
            }

            // increase multidimentional vector with a multidimentional vector
            void operator+=(vector<vector<double>> mat1, const vector<vector<double>>& mat2) {
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1[i].size(); j++) {
                        mat1[i][j] += mat2[i][j];
                    }
                }
            }


            // =============================================                                                                                         
            // Subtraction
            // =============================================  

            // subtraction of a vector and a scalar 
            vector<double> operator-(const vector<double>& vec1, const double& value) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] - value;
                return vec;
            }

            // subtraction of a scalar and a vector 
            vector<double> operator-(const double& value, const vector<double>& vec1) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = value - vec1[i];
                return vec;
            }

            // subtraction of two vectors
            vector<double> operator-(const vector<double>& vec1, const vector<double>& vec2) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] - vec2[i];
                return vec;
            }

            // subtraction of a multidimentional vector and a scalar  
            vector<vector<double>> operator-(const vector<vector<double>>& mat1, const double& value) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = mat1[i][j] - value;
                    }
                }
                return mat;
            }

            // subtraction of a scalar and a multidimentional vector 
            vector<vector<double>> operator-(const double& value, const vector<vector<double>>& mat1) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = value - mat1[i][j];
                    }
                }
                return mat;
            }

            // subtraction of two multidimentional vectors
            vector<vector<double>> operator-(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = mat1[i][j] - mat2[i][j];
                    }
                }
                return mat;
            }

            // decrease vector with a scalar
            void operator-=(vector<double> vec1, const double& value) {
                for (int i{}; i < vec1.size(); i++) vec1[i] -= value;
            }

            // decrease vector with a vector
            void operator-=(vector<double> vec1, const vector<double>& vec2) {
                for (int i{}; i < vec1.size(); i++) vec1[i] -= vec2[i];
            }

            // decrease multidimentional vector with a scalar
            void operator-=(vector<vector<double>> mat1, const double& value) {
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1[i].size(); j++) {
                        mat1[i][j] -= value;
                    }
                }
            }

            // decrease multidimentional vector with a multidimentional vector
            void operator-=(vector<vector<double>> mat1, const vector<vector<double>>& mat2) {
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1[i].size(); j++) {
                        mat1[i][j] -= mat2[i][j];
                    }
                }
            }


            // =============================================                                                                                         
            // Moltiplication
            // =============================================  

            // moltiplication of a vector and a scalar 
            vector<double> operator*(const vector<double>& vec1, const double& value) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] * value;
                return vec;
            }

            // moltiplication of a scalar and a vector 
            vector<double> operator*(const double& value, const vector<double>& vec1) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = value * vec1[i];
                return vec;
            }

            // moltiplication of a vector and a scalar
            void operator*=(vector<double> vec1, const double& value) {
                for (int i{}; i < vec1.size(); i++) vec1[i] *= value;
            }

            // moltiplication of two vectors
            vector<vector<double>> operator*(const vector<double>& vec1, const vector<double>& vec2) {
                vector<vector<double>> mat = zeros(vec1.size(), vec2.size()); 
                for (int i{}; i < vec1.size(); i++) {
                    for (int j{}; j < vec2.size(); j++) {
                        mat[i][j] = vec1[i] * vec2[j];
                    }
                }
                return mat;
            }

            // moltiplication of a multidimentional vector and a scalar
            vector<vector<double>> operator*(const vector<vector<double>>& mat1, const double& value) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = mat1[i][j] * value;
                    }
                }
                return mat;
            }

            // moltiplication of a scalar and a multidimentional vector 
            vector<vector<double>> operator*(const double& value, const vector<vector<double>>& mat1) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = value * mat1[i][j];
                    }
                }
                return mat;
            }

            // moltiplication of a multidimentional vector and a scalar
            void operator*=(vector<vector<double>> mat1, const double& value) {
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat1[i][j] *= value;
                    }
                }
            }

            // moltiplication of two multidimentional vectors
            vector<vector<double>> operator*(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2) {
                vector<vector<double>> mat = zeros(mat1.size(), mat2.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        for (int k{}; k < mat2.front().size(); k++) {
                            mat[i][k] += mat1[i][j] * mat2[j][k];
                        }
                    }
                }
                return mat;
            }

            // moltiplication of a multidimentional vector and a vector
            vector<double> operator*(const vector<vector<double>>& mat1, const vector<double>& vec1) {
                vector<double> vec = zeros(mat1.size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < vec1.size(); j++) {
                        vec[i] += mat1[i][j] * vec1[j];
                    }
                }
                return vec;
            }

            // moltiplication of a vector and a multidimentional vector 
            vector<double> operator*(const vector<double>& vec1, const vector<vector<double>>& mat1) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) {
                    for (int j{}; j < mat1.size(); j++) {
                        vec[i] += vec1[j] * mat1[j][i];
                    }
                }
                return vec;
            }


            // =============================================                                                                                         
            // division
            // =============================================  

            // division of a vector and a scalar 
            vector<double> operator/(const vector<double>& vec1, const double& value) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] / value;
                return vec;
            }

            // division of a scalar and a vector 
            vector<double> operator/(const double& value, const vector<double>& vec1) {
                vector<double> vec = zeros(vec1.size());
                for (int i{}; i < vec1.size(); i++) vec[i] = value / vec1[i];
                return vec;
            }

            // division of a vector and a scalar
            void operator/=(vector<double> vec1, const double& value) {
                for (int i{}; i < vec1.size(); i++) vec1[i] /= value;
            }

            // division of two vectors
            vector<vector<double>> operator/(const vector<double>& vec1, const vector<double>& vec2) {
                vector<vector<double>> mat = zeros(vec1.size(), vec2.size());
                for (int i{}; i < vec1.size(); i++) {
                    for (int j{}; j < vec2.size(); j++) {
                        mat[i][j] = vec1[i] / vec2[j];
                    }
                }
                return mat;
            }

            // division of a multidimentional vector and a scalar
            vector<vector<double>> operator/(const vector<vector<double>>& mat1, const double& value) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = mat1[i][j] / value;
                    }
                }
                return mat;
            }

            // division of a scalar and a multidimentional vector 
            vector<vector<double>> operator/(const double& value, const vector<vector<double>>& mat1) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = value / mat1[i][j];
                    }
                }
                return mat;
            }

            // division of a multidimentional vector and a scalar
            void operator/=(vector<vector<double>> mat1, const double& value) {
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat1[i][j] /= value;
                    }
                }
            }

            // division of two multidimentional vectors
            vector<vector<double>> operator/(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2) {
                vector<vector<double>> mat = zeros(mat1.size(), mat2.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        for (int k{}; k < mat2.front().size(); k++) {
                            mat[i][k] += mat1[i][j] / mat2[j][k];
                        }
                    }
                }
                return mat;
            }

            // division of a multidimentional vector and a vector
            vector<double> operator/(const vector<vector<double>>& mat1, const vector<double>& vec1) {
                vector<double> vec = zeros(mat1.size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < vec1.size(); j++) {
                        vec[i] += mat1[i][j] / vec1[j];
                    }
                }
                return vec;
            }

            // reciprocal of a vector
            vector<double> reciprocal(const vector<double>& vec) {
                vector<double> v = zeros(vec.size());
                for (int i{}; i < v.size(); i++) v[i] = 1. / vec[i];
                return v;
            }

            // reciprocal of a multidimentional vector
            vector<vector<double>> reciprocal(const vector<vector<double>>& mat1) {
                vector<vector<double>> mat = zeros(mat1.size(), mat1.front().size());
                for (int i{}; i < mat1.size(); i++) {
                    for (int j{}; j < mat1.front().size(); j++) {
                        mat[i][j] = 1. / mat1[i][j];
                    }
                }
                return mat;
            }

        } // namespace vectors









        //     // The design requirement is for this to fit in the space of 2 doubles
        //     static_assert(
        //         sizeof(measurement) <= 2 * bitwidth::base_size + sizeof(double),
        //         "Measurement class is too large");

        //     constexpr inline measurement operator*(double val, const units::unit& unit_base) {
        //         return { val, unit_base};
        //     }

        //     constexpr inline measurement operator*(const units::unit& unit_base, double val) {
        //         return { val, unit_base};
        //     }

        //     constexpr inline measurement operator/(double val, const units::unit& unit_base) {
        //         return { val, unit_base.inv() };
        //     }

        //     constexpr inline measurement operator/(const unit& unit_base, double val) {
        //         return { 1.0 / val, unit_base};
        //     }

        // } // namespace measurements
    
    } // namespace tools

} // namespace physics




namespace std {

    // Hash function for unit_data
    template<> struct hash<physics::tools::units::unit_data> {
        size_t operator()(const physics::tools::units::unit_data& x) const noexcept {
            uint32_t val;
            std::memcpy(&val, &x, sizeof(val));
            return hash<uint32_t>()(val);
        }
    };

    // Defining the hash functions for a unit and precise_unit so they can be used
    // in unordered_map
    template<>
    struct hash<physics::tools::units::unit> {
        size_t operator()(const physics::tools::units::unit& x) const {
            return hash<physics::tools::units::unit_data>()(x.base_units()) ^ hash<float>()(x.cround());
        }
    };

    template<>
    struct hash<physics::tools::units::precise_unit> {
        size_t operator()(const physics::tools::units::precise_unit& x) const {
            return hash<physics::tools::units::unit_data>()(x.base_units()) ^ hash<double>()(x.cround());
        }
    };

} // namespace std 


































/*

            class measure {

                public:

                    // =============================================
                    // class members
                    // =============================================
                    
                    double m_value;

                    double m_error;  


                    // =============================================
                    // constructors and destructor
                    // =============================================

                    measure(const double& value, const double& error = 0.) : m_value{value}, m_error{error} {}

                    ~measure() {}


                    // =============================================
                    // set and get methods
                    // =============================================       

                    void set_value(const double& value) { m_value = value; }         

                    void set_error(const double& error) { m_error = error; }         

                    constexpr double get_value() const { return m_value; }
                    
                    constexpr double get_error() const { return m_error; }

                    measure get_measure() const { return *this; }

                    
                    // =============================================
                    // print methods
                    // =============================================   

                    void print() const { std::cout << get_value() << " ± " << get_error(); }   

            };


            class measurement : public measure, public unit {

                public:

                    // =============================================
                    // constructors and destructor
                    // =============================================

                    measurement(const double& value, const double& error, const unit_data& base_unit) : measure(value, error), precise_unit(base_unit) {}
                    
                    // measurement(const double& value, const double& error, const int& __base, const char* __prefix = "") : measure(value, error), unit(__base, __prefix) {}
                    
                    // measurement(const double& value, const double& error, const char* __base, const char* __prefix = "") : measure(value, error), unit(__base, __prefix) {}

                    // measurement(const measurement& m) : measure(m.get_measure()), unit(m.get_unit()) {}

                    ~measurement() {}


                    // =============================================
                    // set and get methods
                    // =============================================       
                    
                    measurement get_measurement() const { return *this; }


                    // // =============================================
                    // // print methods
                    // // =============================================   

                    // void print() const {
                    //     measure::print(); 
                    //     unit::print(); 
                    // }

            };


        } // namespace measurement
        
    } // namespace tools

}

        // namespace measurements {




        //     class measurement : public measure, public unit {

        //         public:

        //             // =============================================
        //             // constructors and destructor
        //             // =============================================

        //             measurement(const double& value, const double& error, const int& __base, const int& __prefix = 6) : measure(value, error), unit(__base, __prefix) {}
                    
        //             measurement(const double& value, const double& error, const int& __base, const char* __prefix = "") : measure(value, error), unit(__base, __prefix) {}
                    
        //             measurement(const double& value, const double& error, const char* __base, const char* __prefix = "") : measure(value, error), unit(__base, __prefix) {}

        //             measurement(const measurement& m) : measure(m.get_measure()), unit(m.get_unit()) {}

        //             ~measurement() {}


        //             // =============================================
        //             // set and get methods
        //             // =============================================       
                    
        //             measurement get_measurement() const { return *this; }


        //             // =============================================
        //             // print methods
        //             // =============================================   

        //             void print() const {
        //                 measure::print(); 
        //                 unit::print(); 
        //             }

        //     };


        // } // namespace measurements



        // namespace constants {

        //     enum class constant_enum { G = 0,  // [m^3 kg^-1 s^-1])
        //                                K = 1 
        //                              };

        //     class constant {
                
        //         protected: 

        //             // =============================================
        //             // convertion methods
        //             // =============================================

        //             constexpr int const_chars_to_int(const char* __prefix) const { 
        //                 if      (__prefix == "G")  { return 0; }
        //                 else if (__prefix == "epsilon_0")  { return 1; }
        //                 else if (__prefix == "K")  { return 2; }
        //                 else {
        //                     std::cerr << "Invalid int for the constant_enum convertion to const char*" << std::endl; 
        //                     exit(-11); 
        //                 }
        //             }

        //             const char* int_to_const_chars(const int& n) const { 
        //                 if      (n == 0)  { return "G"; }
        //                 else if (n == 1)  { return "epsilon_0"; }
        //                 else if (n == 2)  { return "K"; }
        //                 else {
        //                     std::cerr << "Invalid int for the constant_enum convertion to const char*" << std::endl; 
        //                     exit(-11); 
        //                 }
        //             }

        //             constexpr double int_to_double(const int& n) const { 
        //                 if      (n == 0)  { return 6.6743015E-11;      }
        //                 else if (n == 1)  { return 8.854187812813E-12; }
        //                 else if (n == 2)  { return 8.987551792314E19;  } 
        //                 else {
        //                     std::cerr << "Invalid int for the constant_enum convertion to the constant value" << std::endl; 
        //                     exit(-11); 
        //                 }
        //             }

                        
        //         public: 

        //             // =============================================
        //             // class member
        //             // =============================================

        //             constant_enum m_constant; 


        //             // =============================================
        //             // constructors and destructor
        //             // =============================================

        //             constant(const char* __constant) { m_constant = constant_enum(const_chars_to_int(__constant)); }
                                        
        //             constant(const int& __constant) { m_constant = constant_enum(__constant); }


        //             // =============================================
        //             // set and get methods
        //             // =============================================   

        //             inline constexpr constant_enum get_constant() { return m_constant; }

        //             inline constexpr int get_constant_index() const { return static_cast<std::underlying_type<constant_enum>::type>(m_constant); }

        //             inline constexpr double get_constant_value() const { return int_to_double(get_constant_index()); }


        //             // =============================================
        //             // print methods
        //             // =============================================   
                    
        //             inline void print() const { std::cout << "- constant " << int_to_const_chars(get_constant_index()) << " = " << get_constant_value() << std::endl; }

        //     };

        // }


        // namespace position {

        //     using namespace vector_algebra;
        //     using namespace tools::units;

        //     class coordinates : public unit {
                
        //         protected: 

        //             // =============================================
        //             // class members
        //             // =============================================
                
        //             // coordinates:     [x] [y] [z] 
                    
        //             std::vector<double> m_coordinates;


        //         public:  

        //             // =============================================
        //             // constructors and destructor
        //             // =============================================

        //             coordinates() : m_coordinates{zeros(3)}, unit(1, 6) {}

        //             coordinates(const std::vector<double>& coord, const int& __prefix = 6) : m_coordinates{coord}, unit(1, __prefix) {}
                    
        //             coordinates(const std::vector<double>& coord, const char* __prefix = "") : m_coordinates{coord}, unit(1, __prefix) {}

        //             ~coordinates() {}
                    
        //             // =============================================
        //             // set methods
        //             // =============================================

        //             void set_coordinates(const std::vector<double>& coord) { m_coordinates = coord; }
                    
        //             void set_coordinate_x(const double& x) { m_coordinates[0] = x; }

        //             void set_coordinate_y(const double& y) { m_coordinates[1] = y;  }

        //             void set_coordinate_z(const double& z) { m_coordinates[2] = z; }
                
                    
        //             // =============================================
        //             // get methods
        //             // =============================================

        //             // const coordinates get_coordinates() const { return *this; }

        //             std::vector<double> get_coordinates() const { return m_coordinates; }

        //             double get_coordinate_x() const { return m_coordinates[0]; }

        //             double get_coordinate_y() const { return m_coordinates[1]; }

        //             double get_coordinate_z() const { return m_coordinates[2]; }
                    
        //             double get_magnitude() const {
        //                 return sqrt(pow(m_coordinates[0], 2) +                 
        //                             pow(m_coordinates[1], 2) + 
        //                             pow(m_coordinates[2], 2));
        //             }        

        //             double get_distance(const std::vector<double>& coord) const {        
        //                 return sqrt(pow(coord[0] - m_coordinates[0], 2) + 
        //                             pow(coord[1] - m_coordinates[1], 2) + 
        //                             pow(coord[2] - m_coordinates[2], 2)); 
        //             }
                    
        //             double get_rho() const { return sqrt(pow(m_coordinates[0], 2) + pow(m_coordinates[1], 2)); }

        //             double get_phi() const { return atan2(m_coordinates[1], m_coordinates[0]); }     

        //             double get_phi(const std::vector<double>& coord) const { return atan2(coord[1] - m_coordinates[1], coord[0] - m_coordinates[0]); }

        //             double get_theta() const { return acos(m_coordinates[2] / get_magnitude()); }
            
        //             double get_theta(const std::vector<double>& coord) { return acos((coord[2] - m_coordinates[2]) / get_distance(coord)); }

        //             std::vector<double> get_direction() const {
        //                 return { cos(get_phi()), sin(get_phi()), m_coordinates[2] / get_magnitude() };
        //             } 

        //             std::vector<double> get_direction(const std::vector<double>& coord1) const {
        //                 return { cos(get_phi(coord1)), sin(get_phi(coord1)), (coord1[2] - m_coordinates[2]) / get_distance(coord1) };
        //             } 
                    

        //             // =============================================
        //             // print methods
        //             // =============================================

        //             void print() const {
        //                 std::cout << "- coordinates = ";
        //                 for (auto i : m_coordinates) std::cout << "[" << i << "]\t"; 
        //                 unit::print();
        //             }

        //     };

            // class velocity {

            //     protected: 

            //         // =============================================
            //         // class members
            //         // =============================================
                
            //         // velocity:     [x] [y] [z] 
                    
            //         std::vector<double> m_velocity = zeros(3);
              

            //     public:  

            //         // =============================================
            //         // constructors
            //         // =============================================

            //         velocity(const std::vector<double>& vel, const int& __prefix = 6) : m_velocity{vel}, tools::units::unit(1, __prefix) {}
                    
            //         // =============================================
            //         // set methods
            //         // =============================================

            //         void set_velocity(const std::vector<double>& vel) { m_velocity = vel; }

            //         void set_velocity_x(const double& x) { m_velocity[0] = x; }

            //         void set_velocity_y(const double& y) { m_velocity[1] = y;  }

            //         void set_velocity_z(const double& z) { m_velocity[2] = z; }


            //         // =============================================
            //         // get methods
            //         // =============================================

            //         std::vector<double> get_velocity() const { return m_velocity; }

            //         double get_velocity_x() const { return m_velocity[0]; }

            //         double get_velocity_y() const { return m_velocity[1]; }

            //         double get_velocity_z() const { return m_velocity[2]; }
                    
            //         double get_magnitude() const {
            //             return sqrt(pow(m_velocity[0], 2) +                 
            //                         pow(m_velocity[1], 2) + 
            //                         pow(m_velocity[2], 2));
            //         }        

            //         double get_phi() const { return atan2(m_velocity[1], m_velocity[0]); }     

            //         double get_theta() const { return acos(m_velocity[2] / get_magnitude()); }
            
            //         std::vector<double> get_direction() const {
            //             return { cos(get_phi()), sin(get_phi()), cos(get_theta()) };
            //         } 
                

            //         // =============================================
            //         // print methods
            //         // =============================================

            //         void print_velocity() const {
            //             std::cout << "- velocity =    ";
            //             for (auto i : m_velocity) std::cout << "[" << i << "]\t"; 
            //             unit::print();
            //         }

            // };

         // namespace position
     

    


    // namespace objects {

    //     using namespace tools::algebra::vectors;
    //     using namespace tools::measurements;
    //     using namespace tools::position; 


    //     class field :  public coordinates { 

    //         public:

    //             // =============================================
    //             // class members
    //             // =============================================

    //             measurement m_source; 

    //             bool m_field;

    //             std::vector<double> m_attraction;    


    //             // =============================================
    //             // constructor and destructor
    //             // =============================================

    //             field(const measurement& m, const coordinates& coord = coordinates(), const bool& status = false) : m_source{m}, coordinates(coord), m_field{status}, m_attraction{zeros(3)} {}

    //             virtual ~field() {}


    //             // =============================================
    //             // set and get methods
    //             // =============================================

    //             void activate_field() { m_field = true; }

    //             void deactivate_field() { m_field = false; }

    //             void reset_attraction() { m_attraction.clear(); }

    //             void add_attraction(const std::vector<double>& attraction) { m_attraction += attraction; }

    //             std::vector<double> get_attraction() const { return m_attraction; }

    //             virtual std::vector<double> attraction(const coordinates& coord) = 0; 


    //             // =============================================
    //             // print methods
    //             // =============================================   

    //             void print_attraction() const { 
    //                 std::cout << "- attraction = "; 
    //                 for (auto i : get_attraction()) std::cout << "[" << i << "]\t";
    //                 std::cout << std::endl; 
    //             }

    //     };           


    //     class mass : public measurement, public field {
                                
    //         public: 

    //             // =============================================
    //             // constructors and destructor
    //             // =============================================
                                
    //             mass(const double& value, const double& error = 0., const int& __prefix = 6, const coordinates& coord = coordinates(), const bool& field_status = false) : measurement(value, error, 0, __prefix), field(get_measurement(), coord, field_status) {}

    //             mass(const double& value, const double& error = 0., const char* __prefix = "", const coordinates& coord = coordinates(), const bool& field_status = false) : measurement(value, error, 0, __prefix), field(get_measurement(), coord, field_status) {}

    //             ~mass() {}


    //             // =============================================
    //             // print methods
    //             // =============================================
                
    //             void print() const { 
    //                 std::cout << "- mass = "; 
    //                 tools::measurements::measurement::print();
    //             }


    //             // =============================================
    //             // gravitational methods
    //             // =============================================

    //             virtual std::vector<double> attraction(const tools::position::coordinates& coord) override {
    //                 if (m_field == false) {
    //                     std::cout << "Before evaluating the gravitational attraction, you must activate the gravitational field." << std::endl; 
    //                     exit(-11);
    //                 }
    //                 // else if (coord.get_coordinates() == get_coordinates()) { return zeros(3); }
    //                 else return (zeros(3) + 1.); //* (- G * get_value() / pow(get_distance(coord.get_coordinates()), 2)); }
    //             }
            
    //     };

    // } // namespace objects

 // namespace physics
*/
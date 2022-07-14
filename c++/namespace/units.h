
#include <iostream>
#include <cstdint>
#include <cstring>
#include <cmath>


    constexpr int32_t max_neg(uint32_t n_bits) {
        return -(int32_t(1U << (n_bits - 1)));
    }

    namespace bitwidth {
        constexpr uint32_t base_size = sizeof(uint32_t) == 8 ? 8 : 4;
        constexpr uint32_t meter{(base_size == 8) ? 8 : 4};
        constexpr uint32_t second{(base_size == 8) ? 8 : 4};
        constexpr uint32_t kilogram{(base_size == 8) ? 6 : 3};
        constexpr uint32_t ampere{(base_size == 8) ? 6 : 3};
        constexpr uint32_t candela{(base_size == 8) ? 4 : 2};
        constexpr uint32_t kelvin{(base_size == 8) ? 6 : 3};
        constexpr uint32_t mole{(base_size == 8) ? 4 : 2};
    }

    template<typename T>
    constexpr T sqr_power(T a) { return a * a; }

    template<typename X>
    constexpr X power_const_small(X val, int power) {
        return (power == 1) ? val : ((power == -1) ? X(1.0) / val : X(1.0));
    }

    template<typename X>
    constexpr X power_const(X val, int power) {
        return (power > 1) ? sqr_power(power_const(val, power / 2)) * (power % 2 == 0 ? X(1.0) : val) :
            (power < -1) ? X(1.0) / (sqr_power(power_const(val, (-power) / 2)) * ((-power) % 2 == 0 ? X(1.0) : val)) :
                power_const_small(val, power);
    }

    inline float cround(float val) { 
        std::uint32_t bits;
        memcpy(&bits, &val, sizeof(bits));
        bits += 8UL;
        bits &= 0xFFFFFFF0UL;
        memcpy(&val, &bits, sizeof(bits));
        return val;
    }

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

    class unit_data {

        public:

            enum base {
                Meter = 0,
                Second = 1,
                Kilogram = 2,
                Ampere = 3,
                Candela = 4,
                Kelvin = 5,
                Mole = 6,
            };

            static constexpr uint32_t bits[14] = { bitwidth::meter,
                                                   bitwidth::second,
                                                   bitwidth::kilogram,
                                                   bitwidth::ampere,
                                                   bitwidth::candela,
                                                   bitwidth::kelvin,
                                                   bitwidth::mole
                                                };
                

            constexpr unit_data(int meters, 
                                int kilograms, 
                                int seconds, 
                                int amperes, 
                                int kelvins, 
                                int moles, 
                                int candelas) :
                meter_(meters), second_(seconds), kilogram_(kilograms), ampere_(amperes),
                candela_(candelas), kelvin_(kelvins), mole_(moles) {}


            explicit constexpr unit_data(std::nullptr_t) :
                meter_(max_neg(bitwidth::meter)), 
                second_(max_neg(bitwidth::second)),
                kilogram_(max_neg(bitwidth::kilogram)),
                ampere_(max_neg(bitwidth::ampere)),
                candela_(max_neg(bitwidth::candela)),
                kelvin_(max_neg(bitwidth::kelvin)), 
                mole_(max_neg(bitwidth::mole)) {}


            constexpr unit_data operator*(const unit_data& other) const {
                return {
                    meter_ + other.meter_,
                    kilogram_ + other.kilogram_,
                    second_ + other.second_,
                    ampere_ + other.ampere_,
                    kelvin_ + other.kelvin_,
                    mole_ + other.mole_,
                    candela_ + other.candela_,
                };
            }


            constexpr unit_data operator/(const unit_data& other) const {
                return {
                    meter_ - other.meter_,
                    kilogram_ - other.kilogram_,
                    second_ - other.second_,
                    ampere_ - other.ampere_,
                    kelvin_ - other.kelvin_,
                    mole_ - other.mole_,
                    candela_ - other.candela_,
                };
            }


            constexpr unit_data inv() const {
                return {
                    -meter_,
                    -kilogram_,
                    -second_,
                    -ampere_,
                    -kelvin_,
                    -mole_,
                    -candela_};
            }


            constexpr unit_data pow(int power) const { 
                return {
                    meter_ * power,
                    kilogram_ * power,
                    (second_ * power) + rootHertzModifier(power),
                    ampere_ * power,
                    kelvin_ * power,
                    mole_ * power,
                    candela_ * power
                };
            }


            constexpr unit_data root(int power) const {
                return (has_valid_root(power)) ? unit_data( meter_ / power,
                                                          kilogram_ / power,
                                                          second_ / power,
                                                          ampere_ / power,
                                                          kelvin_ / power,
                                                          mole_ / power,
                                                          candela_ / power) : unit_data(nullptr);
            }


            constexpr bool operator==(const unit_data& other) const {
                return equivalent_non_counting(other) && mole_ == other.mole_;         
            }

            constexpr bool operator!=(const unit_data& other) const {
                return !(*this == other);
            }

            constexpr bool has_same_base(const unit_data& other) const {
                return equivalent_non_counting(other) && mole_ == other.mole_;
            }

            constexpr bool equivalent_non_counting(const unit_data& other) const {
                return meter_ == other.meter_ && second_ == other.second_ &&
                    kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                    candela_ == other.candela_ && kelvin_ == other.kelvin_;
            }
            
            constexpr bool empty() const {
                return meter_ == 0 && second_ == 0 && kilogram_ == 0 &&
                    ampere_ == 0 && candela_ == 0 && kelvin_ == 0 && mole_ == 0;
            }

            constexpr int unit_type_count() const {
                return ((meter_ != 0) ? 1 : 0) + ((second_ != 0) ? 1 : 0) +
                    ((kilogram_ != 0) ? 1 : 0) + ((ampere_ != 0) ? 1 : 0) +
                    ((candela_ != 0) ? 1 : 0) + ((kelvin_ != 0) ? 1 : 0) +
                    ((mole_ != 0) ? 1 : 0);
            }
            
            constexpr int meter() const { return meter_; }
            
            constexpr int kg() const { return kilogram_; }
            
            constexpr int second() const { return second_; }
            
            constexpr int ampere() const { return ampere_; }
            
            constexpr int kelvin() const { return kelvin_; }
            
            constexpr int mole() const { return mole_; }
            
            constexpr int candela() const { return candela_; }


        private: 

            constexpr bool has_valid_root(int power) const {
                return meter_ % power == 0 && second_ % power == 0 &&
                    kilogram_ % power == 0 && ampere_ % power == 0 &&
                    candela_ % power == 0 && kelvin_ % power == 0 &&
                    mole_ % power == 0;
            }      

            constexpr int rootHertzModifier(int power) const {
                return (second_ * power == 0 || power % 2 != 0) ? 0 :
                    (power / 2) * ((second_ < 0) || (power < 0) ? 9 : -9);
            }    

            signed int meter_ : bitwidth::meter;
            signed int second_ : bitwidth::second;  
            signed int kilogram_ : bitwidth::kilogram;
            signed int ampere_ : bitwidth::ampere;
            signed int candela_ : bitwidth::candela;  
            signed int kelvin_ : bitwidth::kelvin;
            signed int mole_ : bitwidth::mole;

    };    

    class unit {

        public:
             
            constexpr unit() noexcept {}

            explicit constexpr unit(const unit_data& base_unit) : base_units_(base_unit) {}

            constexpr unit(const unit_data& base_unit, double mult) : base_units_(base_unit), multiplier_(static_cast<float>(mult)) {}

            /// Construct unit from base unit and a multiplier
            constexpr explicit unit(const unit_data& base_unit, float mult) : base_units_(base_unit), multiplier_(mult) {}

            /// Take the double and unit in either order for simplicity
            constexpr unit(double mult, const unit& other) : unit(other.base_units_, mult * other.multiplier()) {}

            /// Unit multiplication
            constexpr unit operator*(const unit& other) const {
                return {base_units_ * other.base_units_, multiplier() * other.multiplier()};
            }

            /// Division operator
            constexpr unit operator/(const unit& other) const {
                return {base_units_ / other.base_units_, multiplier() / other.multiplier()};
            }

            /// Invert the unit (take 1/unit)
            constexpr unit inv() const {
                return {base_units_.inv(), 1.0 / multiplier()};
            }

            /// take a unit to an integral power
            constexpr unit pow(int power) const {
                return unit{base_units_.pow(power), power_const(multiplier_, power)};
            }

            /// Test for unit equivalence to within nominal numerical tolerance (6
            /// decimal digits)
            bool operator==(const unit& other) const {
                if (base_units_ != other.base_units_) { return false; }
                if (multiplier_ == other.multiplier_) { return true; }
                return compare_round_equals(multiplier_, other.multiplier_);
            }

            bool operator!=(const unit& other) const { return !operator == (other); }

            // Test for exact numerical equivalence
            constexpr bool is_exactly_the_same(const unit& other) const {
                return base_units_ == other.base_units_ && multiplier_ == other.multiplier_;
            }

            /// Check if the units have the same base unit 
            constexpr bool has_same_base(const unit& other) const {
                return base_units_.has_same_base(other.base_units_);
            }

            constexpr bool has_same_base(const unit_data& base) const {
                return base_units_.has_same_base(base);
            }

            /// Check if the units have the same base unit  
            constexpr bool equivalent_non_counting(const unit& other) const {
                return base_units_.equivalent_non_counting(other.base_units_);
            }

            constexpr bool equivalent_non_counting(const unit_data& base) const {
                return base_units_.equivalent_non_counting(base);
            }

            /// Check if the units are in some way convertible to one another
            constexpr bool is_convertible(const unit& other) const {
                return base_units_.equivalent_non_counting(other.base_units_);
            }

            /// Check if the base units are in some way directly convertible to one
            /// another
            constexpr bool is_convertible(const unit_data& base) const {
                return base_units_.equivalent_non_counting(base);
            }

            /// Get the number of different base units used
            constexpr int unit_type_count() const {
                return base_units_.unit_type_count();
            }

            /// Extract the base unit Multiplier
            constexpr double multiplier() const {
                return static_cast<double>(multiplier_);
            }

            /// Extract the base unit Multiplier as a single precision float
            constexpr float multiplier_f() const { return multiplier_; }

            constexpr unit_data base_units() const { return base_units_; }


        private:

            friend class precise_unit;

            unit_data base_units_{0, 0, 0, 0, 0, 0, 0};
            
            float multiplier_{1.0};

    };


    constexpr precise_unit hertz(detail::unit_data(0, 0, -1, 0, 0, 0, 0));

    constexpr precise_unit volt(detail::unit_data(2, 1, -3, -1, 0, 0, 0));

    constexpr precise_unit newton(detail::unit_data(1, 1, -2, 0, 0, 0, 0));
    
    constexpr precise_unit pascal(detail::unit_data(-1, 1, -2, 0, 0, 0, 0));
    
    constexpr precise_unit joule(detail::unit_data(2, 1, -2, 0, 0, 0, 0));
    
    constexpr precise_unit watt(detail::unit_data(2, 1, -3, 0, 0, 0, 0));
    
    constexpr precise_unit coulomb(detail::unit_data(0, 0, 1, 1, 0, 0, 0));
    
    constexpr precise_unit farad(detail::unit_data(-2, -1, 4, 2, 0, 0, 0));
    
    constexpr precise_unit ohm(detail::unit_data(2, 1, -3, -2, 0, 0, 0));
        
    constexpr precise_unit weber(detail::unit_data(2, 1, -2, -1, 0, 0, 0));
    
    constexpr precise_unit tesla(detail::unit_data(0, 1, -2, -1, 0, 0, 0));


    constexpr precise_unit Hz = hertz;
    constexpr precise_unit V = volt;
    constexpr precise_unit N = newton;
    constexpr precise_unit Pa = pascal;
    constexpr precise_unit J = joule;
    constexpr precise_unit W = watt;
    constexpr precise_unit C = coulomb;
    constexpr precise_unit F = farad;
    constexpr precise_unit Wb = weber;
    constexpr precise_unit T = tesla;


    // Extra SI units
    constexpr precise_unit bar(100000.0, Pa);
    constexpr precise_unit atm(101325.0, Pa);
    
    // Distance units
    constexpr precise_unit cm(0.01, m);
    constexpr precise_unit km(1000.0, m);
    constexpr precise_unit mm(0.001, m);
    constexpr precise_unit nm(1e-9, m);

    // Volume units
    constexpr precise_unit L{0.001, m* m* m};
    constexpr precise_unit mL{0.001, L};
    // mass units
    constexpr precise_unit g(0.001, kg);
    constexpr precise_unit mg(0.001, g);


    // Time unit
        constexpr precise_unit min(60.0, s);
        constexpr precise_unit ms(0.001, s);
        constexpr precise_unit ns(1e-9, s);
        constexpr precise_unit hr(60.0, min);
        constexpr precise_unit h(60.0, min);
        constexpr precise_unit day(24.0, hr);
        constexpr precise_unit week(7.0, day);
        constexpr precise_unit yr(8760.0, hr);  // median calendar year;
        constexpr precise_unit fortnight(14, day);

        constexpr precise_unit sday{365.24 / 366.24, day};  // sidereal day
        constexpr precise_unit syr(365.256363004, day);  // sidereal year
        constexpr precise_unit at{365.24219, day* eflag};  // mean tropical year
        constexpr precise_unit aj{365.25, day};  // julian year
        constexpr precise_unit ag{365.2425, day};  // gregorian year
        constexpr precise_unit year = yr;  // standard year for SI

        constexpr precise_unit min = min;
    constexpr precise_unit ms = ms;
    constexpr precise_unit ns = ns;
    constexpr precise_unit hr = hr;
    constexpr precise_unit h = h;
    constexpr precise_unit yr = yr;
    constexpr precise_unit day = day;


        constexpr precise_unit deg(constants::pi / 180.0, rad);

        constexpr precise_unit celsius{1.0, K* eflag};
        constexpr precise_unit degC = celsius;

        constexpr precise_unit fahrenheit{5.0 / 9.0, celsius};
        constexpr precise_unit degF = fahrenheit;
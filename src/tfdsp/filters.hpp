/// General linear filtering methods
/// References:
/// [1] : Julius O Smith CCMRA online books https://ccrma.stanford.edu/~jos/
/// [2] ZDF filters : Vadim Zavalishin "The art of VA filter design" rev 2.0 2018
///
#pragma once
#include <cmath>
#include <random>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

/**
 * General linear filtering methods
 * References:
 *     [1] ZDF filters : Vadim Zavalishin "The art of VA filter design" rev 2.0 2018
 *
 */
namespace tfdsp
{
    
/**Prewarping from analog to digital frequency */
template<typename T> inline T AnalogToDigitalFreq_Bilinear(const T samplingFreq, const T fa)
{
    return samplingFreq / PI * std::tan(PI * fa / samplingFreq);
}
template<typename T> inline T DigitalToAnalogFreq_Bilinear(const T samplingFreq, const T fd)
{
    return samplingFreq / PI * std::atan(PI * fd / samplingFreq);
}

template<typename T> class OneSampleDelay
{
     T _s1{}; //Default init to 0 for floating point types

public:
    OneSampleDelay<T>() {}
    inline T operator()(const T x)
    {
        T y = _s1;
        _s1 = x;
        return y;
    }
};
// H(z) = Y/X(z) = (a + z^-1 ) / (1 + a*z^-1)
template<typename T> class FirstOrderAllPass
{
    T _a;
    T _s1{}; //Default init to 0 for floating point types

public:
    FirstOrderAllPass<T>(const T coeff) 
    {
        _a = coeff;
    }
    inline T operator()(const T x)
    {
        T s = x - _a * _s1;
        T y = _a * s + _s1;
        _s1 = s;

        return y;
    }
};

    
 // H(z) = Y/X(z) = ( b0 + b1 * z^-1 + b2 * x z^-2 ) / ( a0 + a1 * z^-1 + a2 * x z^-2 )
// implemented in transposed direct form II ( TDF2 )
template<typename T> class SecondOrderBiquad
{
    std::array<T,3> _a;
    std::array<T,3> _b;
     T _s1{}; //Default init to 0 for numerical types
     T _s2{}; //Default init to 0 for numerical types

public:
    SecondOrderBiquad<T>(const std::array<T,3>& a, const std::array<T,3>& b)
    {
        _a = a;
        _b = b;

        //Normalise coeffs
        _a[0] /= a[0];
        _a[1] /= a[0];
        _a[2] /= a[0];
        _b[0] /= a[0];
        _b[1] /= a[0];
        _b[2] /= a[0];

    }
    inline T operator()(const T x)
    {
        T y =              _b[0] * x + _s1;
        _s1 = -_a[1] * y + _b[1] * x + _s2;
        _s2 = -_a[2] * y + _b[2] * x;
        return y;
    }
    //Critically damped second order butterworth filter
    //Note that fc is normalised frequency between 0 and 1, 1 being the Nyquist limit ( i.e sampling freq / 2 )
    static SecondOrderBiquad<T> ButterworthLowPass(const T fc)
    {
        T c = 1 / std::tan((PI * fc) / 2);
        std::array<T, 3> b{T(1), T(2), T(1)};
        std::array<T, 3> a{T(1 + std::sqrt(2) * c + c * c), T(2 - 2 * c * c), T(1 - std::sqrt(2) * c + c * c)};
        return SecondOrderBiquad(a,b);
    }
};
//First order low pass with zero delay feedback integrator
template<typename T> class FirstOrderLowPassZdf
{
     T _s1{}; //Default init to 0 for floating point types

public:
    FirstOrderLowPassZdf<T>() {}

    //Note that fc is normalised frequency between 0 and 1, 1 being the Nyquist limit ( i.e sampling freq / 2 )
    inline T operator()(const T x, const T fc )
    {
        //g = w^~_c * T /2 = tan( wc *T / 2 ) ( prewarping )
        //i.e g = tan( pi / 2 * f / (fo/2)) = tan(pi / 2 * fc)
        T g = std::tan((PI / 2) * fc);
        T v = (x - _s1) * g / ( 1 + g );
        T y = v + _s1;
        _s1 = y + v;

        return y;
    }
};

//First order high pass with zero delay feedback
template<typename T> class FirstOrderHighPassZdf
{
     T _s1{}; //Default init to 0 for floating point types

public:
    FirstOrderHighPassZdf<T>() {}

    //Note that fc is normalised frequency between 0 and 1, 1 being the Nyquist limit ( i.e sampling freq / 2 )
    inline T operator()(const T x, const T fc )
    {
        //g = w^~_c * T /2 = tan( wc *T / 2 ) ( prewarping )
        //i.e g = tan( pi / 2 * f / (fo/2)) = tan(pi / 2 * fc)
        T g = std::tan((PI / 2) * fc);
        T v = (x - _s1);
        T y = v / (1 + g);
        _s1 = _s1 + y * 2 * g;

        return y;
    }
};

// Non Linear first order low pass with zero delay feedback
// and OTA-style tanh saturation, using its antiderivative with rectangular kernel for antialising,
// and templated number of Newton-Raphson iterations.
 // see [1] 6.13 and [2]
template<unsigned int iterations> class OTAFirstOrderLowPass
{
    //Delayed signals
     double _s1{};
     double _u1{};

public:
    OTAFirstOrderLowPass<iterations>() {}

    //Note that fc is normalised frequency between 0 and 1, 1 being the Nyquist limit ( i.e sampling freq / 2 )
    double operator()(const double x, const double fc)
    {
        //g = w^~_c * T  = 2 * tan( wc *T / 2 ) ( prewarping )
        //i.e g = tan( pi / 2 * f / (fo/2)) = tan(pi / 2 * fc)
        double g = 2. * std::tan((PI / 2.) * fc);

        // y = v + s[1]
        // s = y + v
        // v = g * tanh(u)
        // u = x - y = x - g*tanh(u) - s[1]
        // solve g*tanh(u) + u + -x + s[1] = 0 by Newton Raphson

        //Initial guess is the linearised version
        double u = (x - _s1) / (1 + g);

        for(unsigned int i=0; i < iterations; ++i)
        {
            //u_{n+1}  = u_n - f(u_n) / f'(u_n)
            double tanh_u = std::tanh(u);
            u = u - (u + g*tanh_u -x + _s1) / (1 + g*(1-tanh_u * tanh_u));
        }

        double v = g * std::tanh(u);
        double y = v + _s1;
        _s1 = y + v;

        return y;
    }
};


}

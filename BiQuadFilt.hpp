/***
 * * Filename: BiQuad.hpp
 * *
 * * Description:
 * * This file contains the canonical transposed-direct form II biquad filter.
 * * It includes a helper factory to create specific filter types.
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cmath>
#include <array>
namespace dsp
{
  template<typename T=float>
  class BiQuad
  {
  public: 
    BiQuad(void) noexcept = default; // Default constructor.
    BiQuad(const BiQuad&) noexcept = default; // Copy constructor.
    BiQuad(BiQuad&&) noexcept = default; // Move constructor.
    BiQuad& operator=(const BiQuad&) noexcept = default; // Copy assignment.
    BiQuad& operator=(BiQuad&&) noexcept = default; // Move assignment.
    ~BiQuad() noexcept = default;       // Destructor.
    struct Taps                         // The *Taps of the filter, described by:
    {                                   // The filter structure.
      double b0{0.0f};                   // Feedforward coefficient.
      double b1{0.0f};                   // Feedforward coefficient.
      double b2{0.0f};                   // Feedforward coefficient.
      double a1{0.0f};                   // Feedback coefficient.
      double a2{0.0f};                   // Feedback coefficient.  
    };                                  // Filter *Taps structure.                                
 // Getter and setter for each filter tap.
  //! Set all taps at once
  void Set(const Taps& t) noexcept
  {
    b0 = T(t.b0);
    b1 = T(t.b1);
    b2 = T(t.b2);
    a1 = T(t.a1);
    a2 = T(t.a2);
    z1 = z2 = T(0);
  }

  //! Process one sample
  T ProcessSample(T x) noexcept
  {
    T y = b0 * x + z1;
    z1 = b1 * x + z2 - a1 * y;
    z2 = b2 * x - a2 * y;
    return y;
  }

  //! Return current taps
  Taps GetTaps(void) const noexcept
  {
    return {
      double(b0),
      double(b1),
      double(b2),
      double(a1),
      double(a2)
    };
  }

  //! Convenience setter
  void SetTaps(const Taps& t) noexcept { Set(t); }

  //! Get individual coefficient by index: 0→b0,1→b1,2→b2,3→a1,4→a2
  bool GetCoeff(unsigned i, double& out) const noexcept
  {
    switch(i)
    {
      case 0: out = double(b0); return true;
      case 1: out = double(b1); return true;
      case 2: out = double(b2); return true;
      case 3: out = double(a1); return true;
      case 4: out = double(a2); return true;
      default: return false;
    }
  }

  //! Set individual coefficient by index
  bool SetCoeff(unsigned i, double v) noexcept
  {
    switch(i)
    {
      case 0: b0 = T(v); return true;
      case 1: b1 = T(v); return true;
      case 2: b2 = T(v); return true;
      case 3: a1 = T(v); return true;
      case 4: a2 = T(v); return true;
      default: return false;
    }
  }
  private:
    T b0{0.0f};                            // Feedforward coefficient b0.
    T b1{0.0f};                            // Feedforward coefficient b1.
    T b2{0.0f};                            // Feedforward coefficient b2.
    T a1{0.0f};                            // Feedback coefficient a1.
    T a2{0.0f};                            // Feedback coefficient a2.
    T z1{0.0f};                            // First state variable.
    T z2{0.0f};                            // Second state variable.
  };                                      // Biquad filter class.

// Factory of filter: Butterworth, Chebyshev, Elliptic, and Bessel.
inline BiQuad<float>::Taps ButterworhLP(
  const double fs,                      // The sampling frequency,
  const double fc,                      // Filter cutoff freq.
  const double Q=0.70710678f)           // Butterworth Q factor.
{                                       // ------------- Butterworth ------------- //
  if (fs<=0.0f||fc<=0.0f||fc>=fs*0.5f)  // Is the sampling frequency or cutoff freq invalid?
    return {};                          // Yes, return an empty Taps structure.     
  // -------------------------------- //
  // Pre-Compute the filter characteristics.
  // -------------------------------- //
  /// Pre-warp the frequency.
  double w0=2*M_PI*fc/fs;          // Pre-warp the frequency.
  double cosw0=std::cos(w0);      // Calculate the cosine of w0.
  double sinw0=std::sin(w0);      // Calculate the sine of w0.
  double alpha=sinw0/(2.0*Q);     // Calculate alpha.
  double b0=(1-cosw0)/2.0;        // Calculate b0.
  double b1=1-cosw0;              // Calculate b1.
  double b2=(1-cosw0)/2;          // Calculate b2.
  double a0=1+alpha;              // Calculate a0.
  double a1=-2*cosw0;             // Calculate a1.
  double a2=1-alpha;              // Calculate a2.
  // --------------------------------- //
  // Return the normalized Taps structure.
  // --------------------------------- //
  return {                              // Return the Taps structure.
    .b0=b0/a0,                       // Our b0 coefficient.
    .b1=b1/a0,                       // Our b1 coefficient.
    .b2=b2/a0,                       // Our b2 coefficient.
    .a1=a1/a0,                       // Our a1 coefficient.
    .a2=a2/a0                        // Our a2 coefficient.
  };                                    // Return the Taps structure.
}                                       // ------------- Butterworth ------------- //
  inline BiQuad<float>::Taps ChebyshevLP(
  const double fs,                      // The sampling frequency.
  const double fc,                      // Filter cutoff freq.
  const double ripple=0.1f,             // Chebyshev ripple factor.
  const double Q=0.70710678f)           // Chebyshev Q factor.
{
  if (fs <= 0 || fc <= 0 || fc >= fs*0.5) return {};

  // digital prototype
  double w0=2.0*M_PI*fc/fs;
  double cosw0=std::cos(w0);
  double sinw0=std::sin(w0);
  // Chebyshev epsilon
  double eps=std::sqrt(std::pow(10.0,ripple/10.0) - 1.0);
  // α = sin(w0)·sinh( (1/n)·asinh(1/ε) ), here n=2 → /2
  double v0=std::asinh(1.0/eps) * 0.5;
  double alpha=sinw0*std::sinh(v0);

  // use the Butterworth b-coeffs
  double b0=(1.0-cosw0)*0.5;
  double b1=1.0-cosw0;
  double b2=(1.0-cosw0)*0.5;
  // a-coeffs
  double a0= 1.0+alpha;
  double a1=-2.0*cosw0;
  double a2=1.0-alpha;
  return {
    .b0=b0/a0,
    .b1=b1/a0,
    .b2=b2/a0,
    .a1= a1/a0,
    .a2=a2/a0
  };
}
inline BiQuad<float>::Taps EllipticFirstOrder(
  const double fs,                   // Sampling frequency
  const double fc,                   // Cutoff frequency
  const double /* ripple */ = 0.0    // (unused for 1st order)
)
{
  if (fs <= 0.0 || fc <= 0.0 || fc >= fs * 0.5)
    return {};

  // first-order low-pass prototype H(s)=1/(s+1)
  const double K    = std::tan(M_PI * fc / fs);
  const double norm = 1.0 / (1.0 + K);

  return {
    .b0 =  K * norm,
    .b1 =  K * norm,
    .b2 =  0.0,
    .a1 = (1.0 - K) * norm,
    .a2 =  0.0
  };
}                          // ---------- EllipticFirstOrder ----- //
    // Second order Elliptic Filter - standar biquad (use  Q and fipple values from Elliptic Table.)
    inline BiQuad<float>::Taps Elliptic(
      const double fs,                  // The sampling frequency.
      const double fc,                  // Filter cutoff freq.
      const double ripple_db=0.1f,         // Elliptic ripple factor.
      const double Q=0.70710678f)       // Elliptic Q factor.
      {                                 // ------------- Elliptic ------------- //
        if (fs<=0.0f||fc<=0.0f||fc>=fs*0.5f) // Is the sampling frequency or cutoff freq invalid?
          return {};                    // Yes, return an empty Taps structure.
      // ------------------------------ //
      // Pre-Compute the filter characteristics.
      // ------------------------------ //
      /// Pre-warp the frequency.
        const double K=std::tan(M_PI*fc/fs);// Pre-warp the frequency.
        const double norm=1.0f/(1.0f+K*K);  // Normalize the coefficients.
        const double epsilon=std::sqrt(std::pow(10.0f,ripple_db/10.0f)-1.0f); // Calculate epsilon.
      // ------------------------------ //
      // Return the Taps structure.
      // ------------------------------ //
        return{                         // Return the filter tap structure.
            .b0=(1.0f*epsilon*(K/Q)*K*K)*norm,// Our b0 coefficient.
            .b1=2.0f*(K*K-1.0f)*norm,   // Our b1 coefficient.
            .b2=(1.0f-epsilon*(K/Q)+K*K)*norm,// Our b2 coefficient.
            .a1=2.0f*(K*K-1.0f)*norm,   // Our a1 coefficient.
            .a2=(1.0f-(K/Q)*K*K)*norm   // Our a2 coefficient.
        };                              /// Return the Taps structure.
    }                                   // ------------- Elliptic ------------- //
  
  // Cascade example.
  inline std::array<BiQuad<float>::Taps, 2> EllipticThirdOrder(
    const double fs,                    // The sampling frequency.
    const double fc,                    // Filter cutoff freq.
    const double ripple_db=0.1f,        // Elliptic ripple factor.
    const double Q=0.70710678f)         // Elliptic Q factor.
  {                                     // -------------- EllipticThirdOrder ------------- //
    // -------------------------------- //
    // Forward it to internal building blocks.
    // -------------------------------- //
    return{                             // Return array of filter structures.
      EllipticFirstOrder(fs,fc,ripple_db),// First order filter.
      Elliptic(fs,fc,ripple_db,Q),     // Second order filter.
    };                                 // Return filter bank,
  }                                    // -------------- EllipticThirdOrder ------------- //
  inline BiQuad<float>::Taps BesselFirstOrder(
  const double fs,                      // The sampling frequency.
  const double fc)                      // Filter cutoff freq.
  {                                     // ----------- BesselFirstOrder ----------- //
    if (fs<=0.0f||fc<=0.0f||fc>=fs*0.5f) // Is the sampling frequency or cutoff freq invalid?
      return {};                        // Yes, return an empty Taps structure.
    // -------------------------------- //
    // Pre-Compute the filter characteristics.
    // -------------------------------- //
    /// Pre-warp the frequency.
    const double K=std::tan(M_PI*fc/fs);// Pre-warp the frequency.
    const double norm=1.0f/(1.0f+K);    // Normalize the coefficients.
    // -------------------------------- //
    // Return the Taps structure.
    // -------------------------------- //
    return {                            // Return the Taps structure.
      .b0=K+norm,                       // Our b0 coefficient.
      .b1=K+norm,                       // Our b1 coefficient.
      .b2=0.0f,                         // Our b2 coefficient.
      .a1=(K-1)*norm,                   // Our a1 coefficient.
      .a2=0.0f                          // Our a2 coefficient.
    };                                  // Return first order Bessel block.
  }                                     // ------------- Bessel ------------- //
  inline BiQuad<float>::Taps Bessel(
      const double fs,                  // The sampling frequency.
      const double fc,                  // Filter cutoff freq.
      const double Q=0.70710678f)       // Bessel Q factor.
  {                                     // ----------- BesselFirstOrder ----------- //
    if (fs<=0.0f||fc<=0.0f||fc>=fs*0.5f) // Is the sampling frequency or cutoff freq invalid?
      return {};                        // Yes, return an empty Taps structure.
    // -------------------------------- //
    // Pre-Compute the filter characteristics.
    // -------------------------------- //
    /// Pre-warp the frequency.
    const double K=std::tan(M_PI*fc/fs);
    const double A0=2.324;              // DC gain normalization
    const double a1_a=2.648;            // a1 coefficient for Bessel filter.
    const double a2_a=2.324;            // a2 coefficient for Bessel filter.
    const double D0 =a2_a+a1_a*K+K*K;   // D0 coefficient for Bessel filter.
    // -------------------------------- //
    // Return the Taps structure.
    // -------------------------------- //
    return {
    .b0=A0*(K*K)/D0,
    .b1=2*A0*(K*K-1.0)/D0,
    .b2=A0*(1.0)/D0,
    .a1=(2*(K*K-1.0))/D0,
    .a2=(1-a1_a*K+a2_a)/D0
  };
  }                                     // ----------- Bessel ------------------ //`
  // Third order Bessel filter cascade example function.
  inline std::array<BiQuad<float>::Taps,2> BesselThirdOrder(
    const double fs,                    // The sampling frequency.
    const double fc,                    // Filter cutoff freq.
    const double Q2=0.6910f)            // Canonical Q for 3rd Order Bessel filter.
  {                                     // ------------- BesselThirdOrder ------------- //
    // -------------------------------- //
    // Forward it to internal building blocks.
    // -------------------------------- //
    return{                             // Return array of filter structures.
      BesselFirstOrder(fs,fc),          // First order filter.
      Bessel(fs,fc,Q2)                  // Second order filter.
    };                                  // Return filter bank,
  }                                     // ------------- BesselThirdOrder ------------- //
}

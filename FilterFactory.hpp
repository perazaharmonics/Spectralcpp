  /*
 * * Filename: FilterFactory.hpp
 * *
 * * Description:
 * * This file contains the API to request filters of either
 * * OnePole, ADSR or BiQuad type.
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include "OnePole.hpp"
#include "BiQuad.hpp"
#include "EnvelopeADSR.hpp"
namespace dsp
{
  enum class FilterType
  {
    FirstOrderLowPass,                  // First order lowpass filter.
    FirstOrderHighPass,                 // First order highpass filter.
    ChebyshevLowPass,                   // Chebyshev lowpass filter.
    ButterworthLP,                      // Butterworth lowpass filter.
    ButterworthHP,                      // Butterworth highpass filter.
    BandPass,                           // Bandpass filter.
    Notch,                              // Notch filter.
    LowShelf,                           // Low shelf filter.
    HighShelf,                          // High shelf filter.
    PeakingEQ,                          // Peaking EQ filter.
    EllipticFirstOrder,                 // Elliptic lowpass filter.
    Elliptic,                           // Elliptic highpass filter.
    EllipticThirdOrder,
    BesselFirstOrder,                   // Bessel first order filter.
    Bessel,                             // Bessel filter.
    BesselThirdOrder,                   // Bessel third order filter.
    ADSR,                               // ADSR envelope generator.
  };                                    // Filter types.
  // Factory function to create a filter of the specified type.
  template<typename T=float>
  class FilterFactory
  {
  public:
    // OnePole filters
    static OnePole<T> OnePoleLP(double fs,double fc)
    {                                   // ------- OnePoleLP --------- //
       OnePole<T> h;                    // Create a OnePole filter object.
       h.SetConf(dsp::OnePole<T>::Conf::Lowpass);// Configure filter topology as LP.
       h.Prepare(fs,fc);                // Prepare the filter with the sample rate and cutoff frequency.
       return h;                        // Return the filter object.
    }                                   // ------- OnePoleLP --------- //
    static OnePole<T> OnePoleHP(double fs,double fc)
    {                                   // ------- OnePoleHP --------- //
       OnePole<T> h;                    // Create a OnePole filter object.
       h.SetConf(dsp::OnePole<T>::Conf::Highpass);// Configure filter topology as HP.
       h.Prepare(fs,fc);                // Prepare the filter with the sample rate and cutoff frequency.
       return h;                        // Return the filter object.
    }                                   // ------- OnePoleHP --------- //
    // BiQuad filters
    static dsp::BiQuad<T> ButterWorthLP(double fs,double fc,double Q)
    {                                    // ------- ButterworthLP --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::ButterworthLP(fs,fc,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- ButterworthLP --------- //
    static dsp::BiQuad<T> ButterworthHP(double fs,double fc,double Q)
    {                                    // ------- ButterworthHP --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::ButterworthHP(fs,fc,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- ButterworthHP --------- //
    static dsp::BiQuad<T> ChebyshevLP(double fs,double fc,double ripple=0.1f,double Q=0.70710678f)
    {                                    // ------- ChebyshevLP --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::ChebyshevLP(fs,fc,ripple,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- ChebyshevLP --------- //
    static dsp::BiQuad<T> ChebyshevHP(double fs,double fc,double ripple=0.1f,double Q=0.70710678f)
    {                                    // ------- ChebyshevHP --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::ChebyshevHP(fs,fc,ripple,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- ChebyshevHP --------- //
    static dsp::BiQuad<T> EllipticFirstOrder(double fs,double fc)
    {                                    // ------- EllipticFirstOrder --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::EllipticFirstOrder(fs,fc));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- EllipticFirstOrder --------- //
    static dsp::BiQuad<T> Elliptic(double fs,double fc,double ripple=0.1f,double Q=0.70710678f)
    {                                    // ------- Elliptic --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::Elliptic(fs,fc,ripple,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- Elliptic --------- //
    static std::array<dsp::BiQuad<T>,2> EllipticThirdOrder(double fs,double fc,double ripple=0.1f,double Q=0.70710678f)
    {                                    // ------- EllipticThirdOrder --------- //
      auto taps_array=dsp::EllipticThirdOrder(fs,fc,ripple,Q); // Get the Elliptic taps.
      std::array<dsp::BiQuad<T>,2> h;    // Filter bank for Elliptic filter.
      h[0].SetTaps(taps_array[0]);       // Set the first filter taps.
      h[1].SetTaps(taps_array[1]);       // Set the second filter taps.
      return h;                          // Return the filter bank.
    }                                    // ------- EllipticThirdOrder --------- //
    static dsp::BiQuad<T> BesselFirstOrder(double fs,double fc)
    {                                    // ------- BesselFirstOrder --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::BesselFirstOrder(fs,fc));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- BesselFirstOrder --------- //
    
    static dsp::BiQuad<T> Bessel(double fs,double fc,double Q=0.70710678f)
    {                                    // ------- Bessel --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::Bessel(fs,fc,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- Bessel --------- //
    static std::array<dsp::BiQuad<T>, 2> BesselThirdOrder(double fs,double fc,double Q=0.70710678f)
    {                                    // ------- BesselThirdOrder --------- //
      auto taps_array=dsp::BesselThirdOrder(fs,fc,Q); // Get the Bessel taps.
      std::array<dsp::BiQuad<T>,2> h;    // Filter bank for Bessel filter.
      h[0].SetTaps(taps_array[0]);       // Set the first filter taps.
      h[1].SetTaps(taps_array[1]);       // Set the second filter taps.
      return h;                          // Return the filter bank.
    }                                    // ------- BesselThirdOrder --------- //
    static dsp::BiQuad<T> AllPass(double fs, double fc, double Q=0.70710678f)
    {                                    // ------- AllPass --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::AllPass(fs,fc,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- AllPass --------- //
    static dsp::BiQuad<T> BandPass(double fs, double fc, double Q=0.70710678f)
    {                                    // ------- BandPass --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::BandPass(fs,fc,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- BandPass --------- //
    static dsp::BiQuad<T> Notch(double fs, double fc, double Q=0.70710678f)
    {                                    // ------- Notch --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::Notch(fs,fc,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- Notch --------- //
    static dsp::BiQuad<T> LowShelf(double fs, double fc, double gain, double Q=0.70710678f)
    {                                    // ------- LowShelf --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::LowShelf(fs,fc,gain,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- LowShelf --------- //
    static dsp::BiQuad<T> HighShelf(double fs, double fc, double gain, double Q=0.70710678f)
    {                                    // ------- HighShelf --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::HighShelf(fs,fc,gain,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- HighShelf --------- //
    static dsp::BiQuad<T> PeakingEQ(double fs, double fc, double gain, double Q=0.70710678f)
    {                                    // ------- PeakingEQ --------- //
      dsp::BiQuad<T> h;                       // Create a BiQuad filter object.
      h.SetTaps(dsp::PeakingEQ(fs,fc,gain,Q));// Set the filter taps.
      return h;                          // Return the filter object.
    }                                    // ------- PeakingEQ --------- //
    // ADSR envelope generator
    static dsp::EnvelopeADSR<T> ADSR(double attack, double decay, double sustain, double release, double sampleRate)
    {                                    // ------- ADSR --------- //
      dsp::EnvelopeADSR<T> adsr;              // Create an ADSR envelope generator object.
      adsr.SetAttack(attack);            // Set the attack time.
      adsr.SetDecay(decay);              // Set the decay time.
      adsr.SetSustain(sustain);          // Set the sustain level.
      adsr.SetRelease(release);          // Set the release time.
      adsr.SetSampleRate(sampleRate);    // Set the sample rate.
      return adsr;                       // Return the ADSR object.
    }                                    // ------- ADSR --------- //

  };
} // namespace dsp

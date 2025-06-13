/***
 * * Filename: OnePole.hpp
 * *
 * * Description:
 * * This file contains the one-pole filter class
 * * A one pole filter is a simple lowpass/highpass filter.
 * * It uses thr bilinear transform to get the filter taps.
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cmath>
#include <cstddef>
#include "spectral/FCWTransforms.h"
namespace dsp
{
  using namespace dsp::spectral; // Use the spectral operations namespace.
  template<typename T=float>
  class OnePole
  {
    public:
      OnePole(void) = default; // Copy constructor.
      OnePole(const OnePole&) noexcept = default; // Copy constructor.
      OnePole(OnePole&&) noexcept = default; // Move constructor.
      OnePole& operator=(const OnePole&) noexcept = default; // Copy assignment.
      OnePole& operator=(OnePole&&) noexcept = default; // Move assignment.
      ~OnePole() noexcept = default;    // Destructor.
    enum class Conf { Lowpass, Highpass }; // Filter Conf.
    // Call this once - *After* choosing the filter Conf (LP or HP).
    bool Prepare(                     // Prepare the filter for processing.
      double fs,                        // The sampling frequency.
      double cutoff) noexcept   // The default filter Conf.
      {                                 // ----------- Prepare ----------- //
        if (fs<=0.0||cutoff<=0.0||cutoff>=fs*0.5)// No fs or cutoff?
          return false;                 // Don't know where to operate, return false.
        this->fs=fs;                    // Set the sampling frequency.
        this->cutoff=cutoff;            // Set the cutoff frequency.
        calcCoeffs();                   // Calculate the filter coefficients.
        this->z1=T(0);                  // Initialize the filter state to zero.
        return true;                    // Return true, we are ready.
      }                                 // ----------- Prepare ----------- //
      T ProcessSample(                  // Process a single sample.
        T x)                            // The sample to process.
       noexcept                         // No exceptions.
       {                                // ----------- ProcessSample --------- //
          T y=static_cast<T>(this->a0*x+this->b1*this->z1);// Apply the filter.
          this->z1=y;                   // Update the filter state.
          return y;                     // Return the filtered sample.
       }                                // ----------- ProcessSample --------- //
      void Reset(void) noexcept { this->z1=0.0; } // Reset the filter state.
      // Getters and setters for the filter parameters.
      inline void SetConf(Conf t) noexcept { this->mode=t;calcCoeffs(); } // Set the filter Conf.
    private:
      Conf mode{Conf::Lowpass};         // The innate filter description.
      T a0{1.0};                        // The innate filter gain.
      T b1{0.0};                        // The feedback coefficient.
      T z1{0.0};                        // The innate filter form of being (state).
      double fs{48000.0};               // The sampling frequency.
      double cutoff{1000.0};            // The cutoff frequency in Hz.
      void calcCoeffs(void) noexcept
      {
        const double phase=std::exp(-2.0*M_PI*this->cutoff/this->fs);
        if (this->mode==Conf::Lowpass)  // Does the user want a lowpass?
        {
          this->a0=1.0-phase;           // Set the gain for lowpass.
          this->b1=phase;               // Set the feedback coefficient for lowpass.
        } 
        else                            // Else they want a highpass.
        {
          this->a0=1.0+phase;           // Set the gain for highpass.
          this->b1=phase-1.0;           // Set the feedback coefficient for highpass.
        }                               // Done setting the coeffs.
      }
    };
}

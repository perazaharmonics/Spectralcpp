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
#include "FCWTransforms.h"
namespace dsp
{
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
    enum class Conf { LowPass, Highpass }; // Filter Conf.
      bool Prepare(                     // Prepare the filter for processing.
      double fs,                        // The sampling frequency.
      double cutoff) noexcept   // The default filter Conf.
      {                                 // ----------- Prepare ----------- //
        if (fs<=0.0||cutoff<=0.0||cutoff>=fs*0.5)// No fs or cutoff?
          return false;                 // Don't know where to operate, return false.
        double phase{std::exp(-2.0*M_PI*cutoff/fs)};// Set the phase factor.
        // ---------------------------- //
        // Now we will arrange the filter coeffiecient based on the Conf of filter.
        // ---------------------------- //
        if (this->mode==Conf::LowPass)    // Do we want a lowwpass filter?
        {                               // Set the filter weights...
           this->a0=1.0-phase;          // a0 is the gain.
           this->b1=phase;              // b1 is the feedback coefficient or zero.
        }                               // Done with lowpass filter
        else                            // Else they want w highpass filter.
        {                               // Set the filter coeffs...
          this->a0=(1.0+phase);         // a0 is the gain.
          this->b1=phase-1.0;           // b1 is the feedback coefficient.
        }                               // Done with highpass filter.
        this->z1=0.0;                   // Reset the filter state.
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
      inline void SetConf(Conf t) noexcept { this->mode=t; } // Set the filter Conf.
      inline Conf GetConf(void) const noexcept { return this->mode; } // Get the filter Conf.
      inline T GetGain(void) const noexcept { return this->a0; } // Get the filter gain.
      inline T GetFeedback(void) const noexcept { return this->b1; } // Get the feedback coefficient.
      inline T GetState(void) const noexcept { return this->z1; } // Get the filter state.
      inline void SetState(T z) noexcept { this->z1=z; } // Set the filter state.
      
    private:
      Conf mode{Conf::LowPass};         // The innate filter description.
      T a0{1.0};                        // The innate filter gain.
      T b1{0.0};                        // The feedback coefficient.
      T z1{0.0};                        // The innate filter form of being (state).

    };
}

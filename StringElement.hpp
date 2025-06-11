 /* 
 * * 
 * *
 * * Filename: StringElement.hpp
 * * Description:
 * * This file contains a StringElement class that represent an ideal
 * * or damped string element in a waveguide network.
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include "DelayBranch.hpp"
#include "FilterFactory.hpp"
namespace dsp::wg
{
  struct StringElement final:Node
  {
    // Reinitialize sample rate, fundamental, and cutoff.
    bool Prepare(double fs, double f0,double fc)
    {                                   // ----------- Prepare ---------------- //
      this->fs=fs;
      this->f0=f0;                      // Fundamental frequency.
      this->fc=fc;                      // Cutoff frequency.
      size_t len=static_cast<size_t>(fs/(2*f0));// Half period of the string element.
      if (len<1)len=1;                  // Ensure minimum length of 1 sample.
      fwd.SetDelay(len);                // Set the delay for the forward branch.
      rev.SetDelay(len);                // Set the delay for the reverse branch.
      this->damp.Prepare(fs,fc);        // Prepare the damping filter with the new parameters.
      this->pos=0;                     // Reset the position to 0.
      return true;                      // Return true if prepared successfully.
    }                                   // ----------- Prepare ---------------- //
    // Write excitation (e.g, noise busrt) into both travelling waves
    void Excite(Sample<float> s) noexcept { fwd.Write(s);rev.Write(-s); } // Write excitation to both branches.
    // Propagate the samples through the 2-port string element.
    void Propagate(size_t n) noexcept override
    {                                   // ----------- Propagate ---------------- //
      for (size_t i=0;i<n;++i)          // For each sample to propagate...
      {
        float inFwd=this->rev.Read(); // Fwd travelling wave meets the reverse travelling wave.
        float inRev=this->fwd.Read(); // Rev travelling wave meets the forward travelling wave.
        float outF=this->damp.ProcessSample(-inFwd); // Apply damping to the forward branch.
        float outR=this->damp.ProcessSample(-inRev); // Apply damping to the reverse branch.
        this->fwd.Write(outF);        // Write the processed sample back to the forward branch.
        this->rev.Write(outR);        // Write the processed sample back to the reverse branch.
        this->pos=outF+outR;          // Update the position sample as the sum of both branches.
      }
    }                                   // ----------- Propagate ---------------- //
    // Output the current pick position sample (output).
    inline float Output(void) const noexcept { return pos; } // Get the current position sample.
  
  private:
    DelayBranch<> fwd,rev;              // Forward and reverse branches.
    OnePole<Sample<float>> damp;               // Damping filter for the string element.
    Sample<float> pos=0;                       // Pick position in the string element.
    double fs=44100.0;                  // Sampling frequency.
    double f0=440.0;                    // Fundamental frequency of the string element.
    double fc=1000.0;                   // Cutoff frequency for the damping filter.  
  };
}

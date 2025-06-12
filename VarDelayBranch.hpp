 /* 
 * * Description:
 * * This file contains a delay line class used to represent a sample
 * * of sound travelling through a waveguide in either fwd (S+) or bwd (S-)
 * * direction. Two DelayBranch objects in opposite directions form one
 * * bidirectional tube/string. Fractional for smoother transitions and pitch bends.
 * *
 * * Author:
 * * JEP J. Enrique Peraza, P&G Labs, LLC
 * *
 */
#pragma once
#include "VarFracDelay.hpp"
#include "WGTypes.hpp"
#include "InstrumentGraph.hpp"

namespace dsp::wg
{
  template<size_t MaxLen = 1 << 15>
  struct VarDelayBranch final:public Node
  {
    VarDelayBranch(void) noexcept=default;
    bool SetDelay(size_t d) noexcept { return (this->fracdl.SetDelay(d)); }
    // Read a fractionally interpolated sample from the head.
    float Read(void) const noexcept { return this->fracdl.Read(); }
    // Look at head without advancing the delay line.
    inline T Peek(void) const noexcept { return this->fracdl.PeekHead(); } // Read the sample at the head of the delay line without advancing it.
    // Ramp the *total* delay (integer+fractional) to 'newDelay' over 'n' samples.
    void RampTo(float newDelay, size_t n) noexcept
    {                                   // ----------- RampTo ----------- //
      this->fracdl.RampTo(newDelay,n);  // Ramp the delay line to the new delay over 'n' samples.
    }                                   // ----------- RampTo ----------- //
    // Advance the fractional-delay pointer by one sample
    inline void Tick(void) noexcept { this->fracdl.Tick(); }
    // Write a sample to tail of the delay line:
    inline void Write(wg::Sample<float> s) noexcept { this->fracdl.Write(s); }
    // PStandard Node::Propagate: shift the contents along.
    void Propagate(size_t n) noexcept override
    {                                   // ----------- Propagate ----------- //
        for (size_t i=0;i<n;++i)        // For each sample to propagate...
        {
          //this->fracdl.Write(this->fracdl.Read());// Propagate the sample across the line.
          T v=fracdl.Read();           // Read the sample from the delay line.
          this->fracdl.Write(v);       // Write the sample back to the delay line.
          fracdl.Tick();               // Advance the delay line by one sample.
        }                              //
    }                                   // ----------- Propagate ----------- //
    void Clear(void) noexcept { this->fracdl.Clear(); } // Clear the delay line.
    inline size_t GetDelay(void) const noexcept { return this->fracdl.GetDelay(); }
    inline size_t GetMaxLen(void) const noexcept { return MaxLen; }
    private:
      VarFracDelay<T,MaxLen> fracdl;// The delay line object.
      //static_assert(MaxLen > 0 && (MaxLen & (MaxLen - 1)) == 0, "MaxLen must be a power of two for efficient wrap-around.");
      //double len{MaxLen};             // The length of the delay line in samples.
   };
} // namespace dsp::wg

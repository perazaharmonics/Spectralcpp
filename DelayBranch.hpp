 /* 
 * * Description:
 * * This file contains a delay line class used to represent a sample
 * * of sound travelling through a waveguide in either fwd (S+) or bwd (S-)
 * * direction. Two DelayBranch objects in opposite directions form one
 * * bidirectional tube/string.
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include "DelayLine.hpp"
#include "WGTypes.hpp"

namespace dsp::wg
{
  template<size_t MaxLen = 1 << 15>
  struct DelayBranch final:public Node
  {
    DelayBranch(void) noexcept=default;
    bool SetDelay(size_t d) noexcept { return (this->dl.SetDelay(d)); }
    float Read(void) const noexcept { return this->dl.Read(); }
    // Write a sample to tail of the delay line:
    inline void Write(wg::Sample<float> s) noexcept { this->dl.Write(s); }
    // Propagate the delay line by 'n' samples.
    // repeatedly read the head and write it back to the tail.
    void Propagate(size_t n) noexcept override
    {                                   // ----------- Propagate ----------- //
        for (size_t i=0;i<n;++i)        // For each sample to propagate...
        {
          //this->dl.Write(this->dl.Read());// Propagate the sample across the line.
          float s=this->dl.Read(); // Read the sample from the delay line.
          this->dl.Write(s);       // Write the sample back to the delay line.
        }  
    }                                   // ----------- Propagate ----------- //
    inline size_t GetDelay(void) const noexcept { return this->dl.GetDelay(); }
    inline size_t GetMaxLen(void) const noexcept { return MaxLen; }
    float Peek(void) const noexcept
    {
      return this->dl.PeekTail(); // Read the sample at the head of the delay line without advancing it.
    }
    private:
      dsp::DelayLine<wg::Sample<float>,MaxLen> dl;// The delay line object.
      //static_assert(MaxLen > 0 && (MaxLen & (MaxLen - 1)) == 0, "MaxLen must be a power of two for efficient wrap-around.");
      //double len{MaxLen};             // The length of the delay line in samples.
   };
} // namespace dsp::wg

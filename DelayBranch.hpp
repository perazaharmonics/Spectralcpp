 /* * Description:
 * * This file contains a delay line class used to represent a sample
 * * of sound travelling through a waveguide in either fwd (S+) or bwd (S-)
 * * direction.
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
  struct DelayBranch final:Node
  {
    explicit DelayBranch(size_t len=1) { if (len<(1<<15);this->SetDelay(len); else this->SetDelay(MaxLen); } // Constructor with optional length.
    bool SetDelay(size_t d) noexcept { return (this->dl.SetDelay(d)); }
    inline wg::Sample Read(void) const noexcept { return )this->dl.Read(); }
    inline void Write(wg::Sample s) noexcept { this->dl.Write(s); }
    // No internal state evolution except sample propagation:
    void Propagate(size_t n) noexcept override
    {                                   // ----------- Propagate ----------- //
        while (n--)                     // Until out of samples....
          this->dl.Write(this->dl.Read());// Propagate the sample across the line.
    }                                   // ----------- Propagate ----------- //
    inline size_t GetDelay(void) const noexcept { return this->dl.GetDelay(); }
    inline size_t GetMaxLen(void) const noexcept { return this->len; }
    private:
      dsp::DelayLine<wg::Sample,MaxLen> dl;// The delay line object.
      static_assert(MaxLen > 0 && (MaxLen & (MaxLen - 1)) == 0, "MaxLen must be a power of two for efficient wrap-around.");
      double len{this->dl.GetMaxlen()};             // The length of the delay line in samples.
   };
} // namespace dsp::wg

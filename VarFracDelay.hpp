 /*
 * Filename: VarFracDelay.hpp
 * 
 * * Description:
 * * A power of two circular buffer whose *read* tap may move smoothly
 * * (for vibrato/pitch bend). Linear interpolation is cheap and stable.
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <array>
#include <cstddef>
#include <cassert>

namespace dsp
{

/// A variable, fractional‐delay ring buffer.
/// MaxLen must be a power of two.
template<typename T = float, std::size_t MaxLen = 1024>
class VarFracDelay {
  static_assert((MaxLen & (MaxLen-1))==0, "MaxLen must be power of two");

public:
  VarFracDelay() { Clear(); }

  /// Set the starting (integer) delay
  /// Must be <= MaxLen-1
  void SetDelay(size_t d) noexcept 
  {
    assert(d<MaxLen);
    delay=float(d);
  }

  /// Write one new sample into the ring.
  void Write(const T& x) noexcept 
  {
    buf[widx]=x;
    if ( ++widx==MaxLen) widx=0;
  }

  /// Read the delayed sample with linear interpolation.
  /// call this *after* Write() / Tick()
  T Read(void) const noexcept 
  {
    // integer & fractional parts of the delay
    float d=delay;
    size_t i0=wid  1;             // last‐written sample
    // -------------------------- //
    // wrap i0
    // -------------------------- //
    if (int(i0)<0)i0+=MaxLen;
    size_t idx0 = (i0+MaxLen-size_t(d))&(MaxLen-1);
    float frac=d-std::floor(d);
    size_t idx1=(idx0+MaxLen-1)&(MaxLen-1);
    return buf[idx0]*(1-frac)+buf[idx1]*frac;
  }

  /// Slide the delay line from its current value
  /// to newDelay over nFrames calls to Tick().
  void RampTo(float newDelay, stdsize_t nFrames) noexcept 
  {
    target=newDelay;
    incr=(target-delay)/float(nFrames);
  }

  /// Call once per sample *after* Write(), to update delay.
  void Tick() noexcept {
    delay += incr;
    // keep it in valid bounds [0, MaxLen-1)
    if (delay < 0.0f) delay+=MaxLen;
    else if (delay>=MaxLen)
      delay=std::fmod(delay, float(MaxLen));
  }

  /// Reset the entire buffer to zero and reset pointers.
  void Clear(void) noexcept 
  {
    buf.fill(T{});
    widx=0;
    delay=0.0f;
    target=0.0f;
    incr=0.0f;
  }

private:
  std::array<T,MaxLen> buf;
  std::size_t  widx=0;
  float delay=0.0f;   // current (fractional) delay in samples
  float target=0.0f;   // next target for ramp
  float incr=0.0f;   // per‐sample increment when ramping
};

} // namespace dsp

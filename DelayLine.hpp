/**
 * Filename: DelayLine.hpp
 * 
 * * Description:
 * * This file containis the implementation of a daly line class.
 * * The delay line is implemented as a circular buffer.
 * * It is used to delay the audio signal by a specified number of samples.
 * * Maxelen is the compile time length in sampples of the delay line.
 * Author:
 * JEP J. Enrique Peraza
 *
 */
#pragma once
#include <array>
#include <cstddef>
#include <cassert>
namespace dsp
{
  // ----------------------------------
  // Fixed size fractional delay line (ring-buffer).
  // Maxelen is the compile time length in samples of the delay line.
  // ----------------------------------
  template<typename T,std::size_t Maxlen>
class DelayLine
{
  public:
    // Get DelayLine size in samples.
    inline size_t GetDelay(void) const noexcept { return this->delay; } // Get the current delay in samples.
    
        // Set integer delay (1 <=d <=maxlen-1)
    bool SetDelay( std::size_t d) noexcept                  // The delay in samples.                    // Set the value of the delay line.
    {                                 // ---------------- SetDelay ---------- 
          if (d==0||d>=this->maxlen)      // Is the number of samples out of range
              return false;               // Yes return false.
          this->delay=d;                  // Set the delay.
          return true;                    // Return true.
    }                                 // ---------------- SetDelay ----------
        // Fast integer delay read and write.
    inline void Write(T x) noexcept 
    {
      buf[this->widx]=x;
      this->Advance();
    }
    const inline T Read(void) const noexcept { return buf[(this->widx+this->maxlen-this->delay)%this->maxlen];}
    // Linear interpolated fractional delay read.
    T ReadFrac(                        
      double N)                       // The fractional delay in sample to read.
      const noexcept                  // No exceptions.
    {                                 // ------------ ReadFrac ------------- //
    // ---------------------------- //
    // Assert that the fractional delay is in range.
    // ---------------------------- //
    assert(N>=0.0&&N<static_cast<float>(this->maxlen-1));
    auto idx=static_cast<size_t>(N);// Get the integer part.
    float frac=N-idx;               // Get the fractional part.
    size_t pos0=(this->widx+this->maxlen-idx)&this->mask;// Get the first index.
    size_t pos1={pos0==0?this->maxlen-1:pos0-1}; // Get the second index.
    // ---------------------------- //
    // Interpolate the two samples. //
    // ---------------------------- //
    return (1.f-frac)*buf[pos0]+frac*buf[pos1];// Return the interpolated value.  
  }                                   // ------------ ReadFrac ------------- //
  inline void Clear(void) noexcept { buf.fill(T{});this->widx=0;this->delay; } // Clear the delay line.
  inline const size_t GetMaxlen(void) const noexcept { return this->maxlen; }
  T PeekTail(void) const noexcept 
  {
    size_t id=(widx+maxlen-1)&mask;
    return buf[id]; // Return the last sample in the buffer.
  }
private:
  const size_t maxlen{Maxlen}; // The maximum length of the delay line.
  static constexpr size_t mask=Maxlen-1;// Mask for the circular buffer.
  std::array<T,Maxlen> buf{};           // The circular buffer.
  size_t widx{0};                       // The write index.
  size_t delay{0};                     // The delay in samples.
  inline void Advance(void) noexcept {widx=(widx+1)&mask;}// Advance the write index.
  static_assert((Maxlen & (Maxlen - 1)) == 0, "maxlen must be a power of two for cheap wrap-around.");
};                                      
} // namespace dsp

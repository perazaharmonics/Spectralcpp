 /*
 * *
 * * Filename: LagrangeFracDelayLine.hpp
 * *
 * * Description:
 * * A time-varying fractional delay using Lagrange interpolation.
 * *  instead of the MLTSine window actiong as the low-pass filter,
 * *  we use a Lagrange FIR Fractional Delay Filter (LFD) to interpolate
 * *  between samples
 * *  Coefficient are computed using the Lagrange polynomial: 
 * *   For k = 0 … N
 * *       hk(μ)=∏m=0m≠kNμ−mk−m,0≤μ<1
 * *       hk​(μ)=m=0m=k​∏N​k−mμ−m​,0≤μ<1
 * *  where mu is the fractional part of the requested delay (delay = D + mu), D=[delay]
 * *  and N is the number of taps of the interpolator,k the index of the tap,
 * *  and m is the index of the sample of the incoming signal.
 * * 
 * * Author:
 * *  JEP J.Enrique Peraza,  
 * *
  */
#pragma once
#include <array>
#include <cstddef>
#include <cassert>
#include <cmath>
#include "DelayLine.hpp"
namespace sig::wg
{
  template<typename T=float,size_t MaxLen=1024>
  class DelayLineLagrange
  {
    // Is the Maxlen ODD!!? No, it has to be even right.....?    
    static_assert((MaxLen&MaxLen-1)==0,"MaxLen must be a power of two");
    public:
      void Clear(void) noexcept
      {                                 // ----------- Clear ----------------- //
        buf.fill(T{});                  // Clear the buffer by filling it with zeros.
        widx=0;                         // Reset the write index to zero.
      }                                 // ----------- Clear ----------------- //
      void Write(const T&x) noexcept
      {                                 // ----------- Write ----------------- //
        this->buf[widx]=x;              // Write the sample x into the buffer at the current write index.
        widx=(widx+1)&(MaxLen-1);       // Increment the write index and wrap it around if it exceeds MaxLen.
      }                                 // ----------- Write ----------------- //
      const T& Peek(size_t idx) const noexcept
      {                                 // ----------- Peek ------------------ //
        // Wrap around if idx is larger than MaxLen                            //
        return buf[(widx-idx)&(MaxLen-1)];// Peek at the sample at index idx from the current write index
      }                                 // ----------- Peek ------------------ //
    
     private:
      std::array<T,MaxLen> buf{};        // Buffer to hold the samples.
      size_t widx=0;                    // Write index for the buffer.
  };
  /// ----------------- Optimal Lagrange fractional-delay FIR ----------------- //
  template<typename T=float,size_t MaxOrder=1024>
  class LagrangeDelayLine
  {
    public:
    void SetOrder(size_t N) noexcept
    {                                   // ----------- SetOrder -------------- //
      order=std::min<size_t>(N,MaxOrder);// Set the order of the Lagrange filter, ensuring it does not exceed MaxLen-1.
    }                                   // ----------- SetOrder -------------- //
    void SetMu(T mu) noexcept
    {                                   // ------------ SetMu ----------------- //
      this->mu=mu-std::floor(mu);       // Set the fractional part of the delay
      build();                          // Rebuild the filter coefficients 
    }                                   // ------------ SetMu ----------------- //
    /// Interpolate sample sitting *id* samples behind the write-ptr of the delay line.
    template<size_t Len>
    T Process(const DelayLine<T, Len>& dl,size_t id) const noexcept
    {                                   // ----------- Process ----------------- //
      T y{};                            // unrolled for small N.
      for (size_t k=0;k<=order;++k)     // For the interpolation order....
        y+=h[k]*dl.Peek(id+k);          // Compute the output sample by summing the products of the coefficients and the corresponding samples from the delay line.
      return y;                         // Return the computed output sample.
    }                                   // ----------- Process ----------------- //
    // Expose the coefficients of the Lagrange filter.
    const std::array<T,MaxOrder+1>& Coeff(void) const noexcept
    {                                   // ----------- Coeff ----------------- //
      return h;                         // Return the coefficients of the Lagrange filter.
    }                                   // ----------- Coeff ----------------- //
    // Expose the order of the Lagrange filter.
    size_t Order(void) const noexcept
    {                                   // ----------- Order ----------------- //
      return order;                     // Return the order of the Lagrange filter.
    }                                   // ----------- Order ----------------- //
    // Expose the fractional part of the delay.
    T Mu(void) const noexcept
    {                                   // ----------- Mu ----------------- //
      return mu;                        // Return the fractional part of the delay.
    }                                   // ----------- Mu ----------------- //
    private:
      void build(void) noexcept
      {                                 // Reconfigure filter topology.
        // For each tap compute h[k]=Π_{m≠k} (μ-m)/(k-m)
        for (size_t k=0;k<=order;++k)   // For each filter coeff to the order of Lagrange...
        {                               // Compute the coefficients for the order
          T c=T(1);                     // Initialize the coefficient to 1.
          for (size_t m=0;m<=order;++m) // For each sample in the filter...
          {                             // Compute the product for the coefficient.
            if (m!=k)                   // If m is not equal to k...
            {                           // Compute the product for the coefficient.
              c*=((mu-T(m))/(T(k)-T(m))); // Compute the product for the coefficient.
            }                           // Done with the product.
          }                             // Done with the coefficient.
          h[k]=c;                       // Store the coefficient in the array.
        }                               // Done with the coefficients.
      }                                 // ----------- build ----------------- //
      
      size_t order=3;                   // Order of the Lagrange filter, default is 3.
      T mu{0};                          // Fractional part of the delay, default is 0.
      std::array<T,MaxOrder+1> h{};     // Coefficients for the Lagrange filter, size is order+1.
  };                                    // ----------- LagrangeDelayLine ----------------- //
  // ----------------- Variable Fractional Delay Line ------------------------------------- //
  template<typename 
    T=float,
    size_t MaxLen=1024,
    size_t MaxOrder=5>
  class LagrangeVarDelayLine
  {
    void Clear(void) noexcept
    {                                   // ----------- Clear ----------------- //
      dl.Clear();                       // Clear the delay line buffer.
      lagrange.SetMu(0);                // Set the fractional part of the delay line to 0.
    }                                   // ----------- Clear ----------------- //
    void Write(const T& x) noexcept
    {                                   // ----------- Write ----------------- //
      dl.Write(x);                      // Write the sample x into the delay line.
    }                                   // ----------- Write ----------------- //
    // Set new delay (can be non-integer) - call at audio-rate or per block.
    void SetDelay(
      T delay,                          // The desired delay in samples.
      size_t order=3) noexcept          // The order of the Lagrange filter
    {                                   // ----------- SetDelay --------------- //
      // Delay greater than, at least 0 and less than Maxlen-1?
      assert(delay>=0&&del<MaxLen-1);     // Ensure delay is bounded within our compact support.
      this->del=delay;                // Set the desired delay.
      lagrange.SetOrder(order);         // Set the order of the Lagrange filter.
      lagrange.SetMu(delay-std::floor(delay)); // Set the fractional part of the delay.
    }                                   // ----------- SetDelay --------------- //
    // Read with optimal Lagrange Fractional Delay Filter (LFD)  interpolator.
    T Read(void) const noexcept
    {                                   // ----------- Read ------------------ //
      const size_t d=static_cast<size_t>(std::floor(del));
      return lagrange.Process(dl,d);   // Read the sample from the delay line using the Lagrange filter.
    }                                   // ----------- Read ------------------ //
    // Convenience smooth ramp in *k* samples.
    void RampTo(                        // Reaches a target delay in *k* samples.
      T d,                              // The target delay in samples.
      size_t k) noexcept                // The number of samples to ramp to the target.
    {                                   // ----------- RampTo ------------------------- //
      incr=(d-del)/static_cast<T>(k); // Compute the increment for the ramp.
      targ=(d);                         // Set the target delay.
    }                                   // ----------- RampTo ------------------------- //
    void Tick(void) noexcept            // Something happened... update myself.
    {                                   // ----------- Tick ------------------ //
      if (incr>0&&del<targ||(incr<0&&del>targ))// We updated the incr or delay greater than target?
      {                                 // Yes, so let's update ourselves....
        del+=incr;                      // We incremented the delay by this much.
        lagrange.SetMu(del-std::floor(del));// So update the interpolator with the delay state.
      }                                 // Done with the update.
    }                                   // ----------- Tick ------------------ //
    // Expose Coefficients
    const std::array<T,MaxOrder+1>& Coeff(
      size_t k) const noexcept          // The index into the coefficients array of the Lagrange filter.
    {                                   // ----------- Coeff ----------------- //
      assert(k<=MaxOrder);              // Ensure the index is within bounds.
      for (size_t i=0;i<k;i++)          // Find index k of the coefficients 
      {                                 // ... of the Lagrange Fractional FIR filter.
        assert(i<MaxOrder+1);           // Ensure the index is within bounds.
        return lagrange.h[i];           // Return the coefficient at index i.
      }                                 // Done getting the coefficient.
    }                                   // ----------- Coeff ----------------- //
    // Expose Process
    T Process(size_t id) const noexcept // Process the delay line with the Lagrange
    {                                   // ----------- Process ----------------- //
      assert(id<MaxLen);                // Ensure the index is within bounds.
      return lagrange.Process(dl,id);   // Process the delay line using the Lagrange
    }                                   // ----------- Process ----------------- //
    // Expose the current delay.
    T Delay(void) const noexcept        // Get the current delay in samples.
    {                                   // ----------- Delay ----------------- //
      return del;                       // Return the current delay in samples.
    }                                   // ----------- Delay ----------------- //
    // Expose the target delay.
    T Target(void) const noexcept       // Get the target delay in samples.
    {                                   // ----------- Target ----------------- //
      return targ;                      // Return the target delay in samples.
    }                                   // ----------- Target ----------------- //
    // Expose the increment for the ramp.
    T Increment(void) const noexcept    // Get the increment for the ramp.
    {                                   // ----------- Increment --------------- //
      return incr;                      // Return the increment for the ramp.
    }                                   // ----------- Increment --------------- //
    private:
      DelayLine<T,MaxLen> dl;           // Delay line to hold the samples.
      LagrangeDelayLine<T,MaxOrder> lagrange; // Lagrange filter
      T del{0};                         // Current delay in samples.
      T targ{0};                        // Target delay in samples.
      T incr{0};                        // Increment for the ramp.    
  };                                    // ----------- LagrangeVarDelayLine ----------------- //
}                                       // namespace sig::wg

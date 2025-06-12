/***
 * * Filename: BodyConvolver.hpp
 * *
 * * Description:
 * * This file contains a partitioned overlap-add FIR
 * * (static IR lightweight convolution static IR <= 4096 samples.)
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <vector>
#include "AudioBuffer.hpp"
namespace dsp::wg
{
  class BodyConvolver
  {
  public:
    // Load an impulse respone both channels must be same length.
    bool LoadImpulseResponse(
      const std::vector<float>& L,
      const std::vector<float>& R)
    {                                   // ----------- LoadImpulseResponse ---- //
      // Is it symmetric, non-empty and less than 4096 samples?
      if (L.size()!=R.size()||L.empty()||L.size()>4096) return false;
      this->irL=L;                      // Load the left channel impulse response.
      this->irR=R;                      // Load the right channel impulse response.
      const size_t N=this->irL.size();  // Get the size of the impulse response.
      this->circ.assign(2*N,0.0f);      // Assign this much to buffer.
      this->widx=0;                     // Reset the write index to zero.
      return true;                      // Always return true, as we assume the IR is valid.
    }                                   // ----------- LoadImpulseResponse ---- //
    // Process one input sample and accumulate into outL/outR.
    void Process(float in, float &outL, float &outR) noexcept
    {                                   // ----------- Process ---------------- //
      const size_t N=irL.size();        // Get the size of the impulse response.
      circ[widx]=in;                    // Write the input sample to the circular buffer at the write index.
      float sumL=0.0f;                  // Initialize the sum for the left channel.
      float sumR=0.0f;                  // Initialize the sum for the right channel.
      size_t idx=widx;                  // Start from the current write index.
      // ------------------------------ //
      // Perform OLA convolution: Walk backwards through the circular buffer
      // and multiply each sample with the corresponding impulse response sample.
      // ------------------------------ //
      for (size_t i=0;i<N;++i)          // For each sample in the impulse response vector...
      {                                 // Perform the convolution.
        sumL+=circ[idx]*irL[i];         // Multiply the circular buffer sample with the left channel impulse response sample and accumulate.
        sumR+=circ[idx]*irR[i];         // Multiply the circular buffer sample with the right channel impulse response sample and accumulate.
        idx=idx?(2*N-1):idx-1;          // Decrement the idx, or set to tail.
      }                                 // Done with the convolution loop.
      // ------------------------------ //
      // Mix into the output channels.
      // ------------------------------ //
      outL+=sumL;                       // Write the accumulated sum to the left output channel.
      outR+=sumR;                       // Write the accumulated sum to the right output channel.
      // ------------------------------ //
      // Advance and wrap write ptr.
      // ------------------------------ //
      if (++widx>=2*N) widx=0;          // 
    }                                   // ----------- Process ---------------- //
    void Reset(void) noexcept           // Reset the circular buffer and write index.
    {                                   // ----------- Reset ------------------ //
      std::fill(circ.begin(), circ.end(), 0.0f); // Clear the circular buffer.
      widx=0;                           // Reset the write index to zero.
    }                                   // ----------- Reset ------------------ //
    // Getters
    size_t GetIRSize(void) const noexcept { return irL.size(); }
    size_t GetBufSize(void) const noexcept { return circ.size(); } // Get the size of the circular buffer.
    size_t GetWriteIndex(void) const noexcept { return widx; }

   private:
     std::vector<float> irL{0.0f};
     std::vector<float> irR{0.0f};         // Impulse response for left and right channels.
     std::vector<float> circ{0.0f};        // Circular buffer for overlap-add convolution.
     size_t widx=0;                        // Write index for the circular buffer.
  };    
}
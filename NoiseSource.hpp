/***
 * * Filename: WhiteNoise.hpp
 * *
 * * Description:
 * * White-Noise source used as the exciter for many systems.
 * *
 * * Author:
 * * JEP J. Enrique Peraza, P&G Labs, LLC
 * *
 */
#pragma once
#include <cstdint>
namespace dsp
{
    // xor-shift white-noise generator
    class WhiteNoise
    {                                   // Class: WhiteNoise
      public:
        explicit WhiteNoise(std::uint32_t seed=1): state(seed) {};
        float Process(void) noexcept
        {                               // ---------- Process ------------- //
          uint32_t x=this->state;       // Get the current state.
          x^=x<<23;x^=x>>17;x^=x<<5;    // Apply the xor-shift algorithm.
          this->state=x;                // Update the state.
          return static_cast<float>(x)*2.3283064365386963e-10f*2.f-1.f; // Scale to [-1,1] range.
        }                               // ---------- Process ------------- //
        inline void Reset(void) noexcept { this->state=1; } // Reset the state to 1.
      private:
        uint32_t state{1};              // The internal state of the generator.
    };                                  // Class: WhiteNoise                                   
}                                       // namespace dsp
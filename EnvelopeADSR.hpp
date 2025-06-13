/***
 * * Filename: EnvelopeADSR.hpp
 * *
 * * Description:
 * * Classic Attack-Decay-Sustain-Release (ASDR) envelope generator.
 * * The time constsants are in seconds. Return the value per-sample
 * * (call inside inner computation loop).
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cstdint>
#include <cmath>
namespace dsp{
  template<typename T=float>
  class EnvelopeADSR
  {
  public:
    bool SetSampleRate(double fs) noexcept {this->fs=fs;return fs>0.0f;}// Set the sample rate in Hz.
    bool SetAttack(double attack) noexcept {this->attack=attack;return attack>0.0f;}// Set the attack time in seconds.
    bool SetDecay(double decay) noexcept {this->decay=decay;return decay>0.0f;}// Set the decay time in seconds.
    bool SetSustain(double sustain) noexcept {this->sustain=sustain;return sustain>=0.0f&&sustain<=1.0f;}// Set the sustain level (0-1).
    bool SetRelease(double release) noexcept {this->release=release;return release>0.0f;}// Set the release time in seconds.
    // Reset the envelope to the initial state.
    void Reset(void) noexcept
    {
      this->state=Stage::Idle;          // Set the state to idle.
      this->val=0.0f;                 // Reset the envelope value.
      this->phase=0.0f;                 // Reset the phase.
    }                                   // Reset the envelope.
    void SetADSR(
        double att,               // The attack time in seconds.
        double dec,               // The decay time in seconds.
        double sus,               // The sustain level (0-1).
        double rel) noexcept      // The release time in seconds.
    {                                   // ------------ SetASDR ------------- //
        this->SetAttack(att);           // Set the attack time.
        this->SetDecay(dec);            // Set the decay time.
        this->SetSustain(sus);          // Set the sustain level.
        this->SetRelease(rel);          // Set the release time.
    }
    void NoteOn(void) noexcept { this->state=Stage::Attack;n=0; } // Start the envelope in attack stage.
    void NoteOff(void) noexcept { this->state=Stage::Release;n=0; } // Start the envelope in release stage.
    // Process the ASDR envelope on this sample.
    T Process(void) noexcept
    {
      switch(this->state)               // According to the state we are in....
      {                                 // Process the samples.
        case Stage::Attack:
          this->val+=inc(this->attack); // Increment the value by the attack increment.
          if (this->val>=1.0f)          // Is the value greater than or equal to 1.0f?
          {                             // Yes, we reached the end of the attack stage.
            this->val=1;                // Clamp the value to 1.
            Next(Stage::Decay);         // Move to the decay stage.
          }                             // Done with attack stage.
          break;                        // Done handling attack, switch to next case.
        case Stage::Decay:              // We are at the decay stage of the machine:
          this->val-=inc(this->decay)*(1.0f*this->sustain); // Decrement the value by the decay increment.
          if (this->val<=this->sustain) // Is the value less than or equal to the sustain level?
          {                             // Yes, we reached the end of the decay stage.
            this->val=this->sustain;    // Clamp the value to the sustain level.
            Next(Stage::Sustain);       // Move to the sustain stage.
          }                             // Done with decay stage.
          break;                        // Done handling decay, switch to next case.
        case Stage::Sustain:            // We are at the sustain stage of the machine:
          val-=inc(this->release)*this->sustain; // Decrement the value by the release increment.
          // In sustain we just hold the value at the sustain level while the note is down.
          // NO-OP, just break.
          break;
        case Stage::Release:            // We are in the release stage
          this->val=inc(this->release);// Graceful tail.
          if (this->val<=0.0f) { this->val=0.0f;Next(Stage::Idle); }
        case Stage::Idle:               // Idla part of the machine:
          default:break;                // Yes, so NO-OP.
        }                               // Done with switch.
        this->n++;                      // Increment the sample counter.
        return static_cast<T>(this->val); // Return the current value of the envelope.
    }                                   // ------------ Process ------------- //
  inline bool IsActive(void) const noexcept { return this->state!=Stage::Idle; }
  inline double GetSampleRate(void) const noexcept { return this->fs; } // Get the sample rate.
  inline double GetAttack(void) const noexcept { return this->attack; } // Get the attack time.
  inline double GetDecay(void) const noexcept { return this->decay; } // Get the decay time.
  inline double GetSustain(void) const noexcept { return this->sustain; } // Get the sustain level.
  inline double GetRelease(void) const noexcept { return this->release; } // Get the release time.

private:
  // Enums to the describe possible state of the envelope state machine.
  enum class Stage {Idle,Attack,Decay,Sustain,Release};
  Stage state{Stage::Idle};             // The innate current state.
  double fs{48000.0};                   // The sample rate in Hz.
  double phase{0.0};                  // The current phase of the envelope.
  // The envelope parameters.
  // The attack, decay, sustain, and release times in seconds.
  double attack{0.01};                  // The attack time in seconds.
  double decay{0.1};                    // The decay time in seconds.
  double sustain{0.7};                  // The sustain level (0-1).
  double release{0.1};                  // The release time in seconds.
  uint64_t n{0};                        // The sample counter.
  T val{0};                             // The current value of the envelope.

  inline constexpr T inc(double time) const noexcept { return (time<=0.0f)?T(1):T(1.0/(time*this->fs));}
  inline void Next(Stage n) noexcept 
  { 
    this->state=n;                      // Set the next state.
    this->n=0;                          // Reset the sample counter.
  }
};
}

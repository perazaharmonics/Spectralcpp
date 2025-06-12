  /*
 * * Filename: BlockBodyConvolver.hpp
 * *
 * * Description:
 * * This file contains a partitioned overlap-add FIR
 * * to process blocks of audio samples. It uses FFT based
 * * convolution to go from O(N) per samples to O(MLogM) per
 * * block.
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>
#include "AudioBuffer.hpp"
#include "FCWTransforms.hpp"
#include "DSPWindows.hpp"
namespace dsp::wg
{
  using namespace dsp::spectral; // Use the spectral operations namespace.
  class BodyConvolverFFT
  {
  public:
    BodyConvolverFFT(void)=default;
    
    // Lifecycle configuration.
    /// Prepare he sample-reate and block size
    /// User must call before LoadImpulseResponse.
    void Prepare(double fs,size_t blockSiz) noexcept
    {                                   // ----------- Prepare ------------------ //
      this->fs=fs;                      // Set the sample rate.
      this->M=blockSiz;                 // Set the block size.
      this->wSiz=int(2*M);              // Set the window size to 2*M.
      this->win.GenerateWindow(Window<float>::WindowType::Hanning,wSiz);// Generate a Hanning window of size 2*M.
      this->engine.SetSampleRate(fs);   // Set the sample rate for the spectral operations engine.
      this->engine.SetSamples(fftsiz);  // Set size of FFT buffers.
    }
    
    // Load a stereo IR (left and right must be same length, >0) and choose
    // a block size M.
    bool LoadImpulseResponse(
        const std::vector<float>& L,
        const std::vector<float>& R)            
    {                                   // -------- LoadImpulseResponse ---- //
        if (L.size()!=R.size()||L.empty()||L.size()>4096) return false; // Check IR validity.
        // ---------------------------- //
        // Copy IR into complex vectors and zero-pad to the next power of 2.
        // ---------------------------- //
        this->irL=L;                    // Load the left channel impulse response.
        this->irR=R;                    // Load the right channel impulse response.
        size_t IRn=L.size();            // Get the size of the impulse response.
        fftsiz=nextPowerOfTwo(std::max<size_t>(wSiz,IRn));
        // Number of M-sample partitions
        part=(IRn+M-1)/M;               // Calculate the number of partitions needed.
        // Allocate and fill partitioned FFTs
        H_L.assign(part,std::vector<std::complex<float>>(fftsiz));
        H_R.assign(part,std::vector<std::complex<float>>(fftsiz));
        std::vector<std::complex<float>> buf(fftsiz);
        // For each partition p, zero-pad M samples into an fftsiz buffer.
        for (size_t p=0;p<part;++p)     // For each partition...
        {
          // left channel.
          std::fill(buf.begin(),buf.end(),std::complex<float>(0,0)); // Zero the buffer for the left channel.
          for (size_t n=0;n<M;++n)
          {
            size_t idx=p*M+n;           // Calculate the index in the impulse response.
            if (idx<IRn)                // Is the index within range?
              buf[n]={L[idx],0.0f};     // Copy the left channel sample into the buffer.
          }                             // Done copying left channel samples.
          H_L[p]=engine.FFT(buf);       // Perform the FFT on the left channel buffer.
          // Right channel.
          std::fill(buf.begin(),buf.end(),std::complex<float>(0,0));// Zero the buffer for the right channel.
          for (size_t n=0;n<M;++n)     // For each sample in the block...
          {                            // Copy the right channel sample into the buffer.
            size_t idx=p*M+n;          // Calculate the index in the impulse response.
            if (idx<IRn)               // Is the index within range?
              buf[n]={R[idx],0.0f};    // Copy the right channel sample into the buffer.
          }                            // Done copying right channel samples.
          H_R[p]=engine.FFT(buf);      // Perform the FFT on the right channel buffer.
        }                              // Done with all partitions.
        // Set the engine to use the FFT size.
        engine.SetSampleRate(fs); // Set the sample rate for the spectral operations engine.
        engine.SetSamples(fftsiz); // Set the size of the FFT buffers.
        // ---------------------------- //
        // Initialize the input buf and accumulators.
        // ---------------------------- //
        inBuf.assign(fftsiz,{0.0f,0.0f}); // Initialize the input buffer with complex zeros.
        accL.assign(fftsiz,{0.0f,0.0f}); // Initialize the left channel accumulator with complex zeros.
        accR.assign(fftsiz,{0.0f,0.0f}); // Initialize the right channel accumulator with complex zeros.
        pidx=0;                        // Reset the partition index to zero.
        return true;                   // Return true if the IR was loaded successfully.
    }                                   // -------- LoadImpulseResponse ---- //
    // Process a block of audio samples.
    // ----------------------------------
    // in: pointer to M mono input samples.
    // outL: pointer to M left output samples.
    // outR: pointer to M right output samples.
    // ----------------------------------
    void ProcessBlock(const float* in, float* outL, float* outR) noexcept
    {                                   // --------- ProcessBlock ------------- //
      assert(M>0&&fftsiz>=2*M);         // Ensure M is set and fftsiz is at least 2*M.
      // ------------------------------ //
      // Build the complex input buffer: copy M samples and window then
      // zero the rest.
      // ------------------------------ //
      for (size_t n=0;n<M;++n)          // For each sample ....
      {                                 // Copy each complex sample.
        inBuf[n]=std::complex<float>(in[n]*win[n],0.0f);
      }                                 // Done copying input samples.
      // Zero-pad the remainder.
      for (size_t n=wSiz;n<fftsiz;++n)  // From the length of the window to fftsiz...
        inBuf[n]={0,0};                     // Zero pad it.
      // ------------------------------ //
      // FFT the input buffer.
      // ------------------------------ //
      auto X=engine.FFT(inBuf);        // Perform the FFT on the input buffer.
      // ------------------------------ //
      // For each IR partition, multiply, IFFT and overlap-add
      // ------------------------------ //
      for (size_t p=0;p<part;++p)          // For each partition...
      {                                 // Overlap and add...
        // ---------------------------- //
        // Point-wise product in frequency domain.
        // ---------------------------- //
        auto YL=engine.PointwiseMul(X,H_L[p]); // Multiply the input spectrum with the left channel partition.
        auto YR=engine.PointwiseMul(X,H_R[p]); // Multiply the input spectrum with the right channel partition.
        // ---------------------------- //
        // Back to time-domain with IFFT.
        // ---------------------------- //
        auto yL=engine.IFFT(YL);        // Perform the IFFT on the left channel partition.
        auto yR=engine.IFFT(YR);        // Perform the IFFT on the right channel partition.
        // ---------------------------- //
        // Overlap-add into accumulators.
        // ---------------------------- //
        for (size_t n=0;n<fftsiz;++n)   // For each segment in the fft sized signal...
        {                               // Overlap-add the result.
          accL[n]+=yL[n];               // Add the left channel segment to the accumulator.
          accR[n]+=yR[n];               // Add the right channel segment to the accumulator.
        }                               // Done with the overlap-add.
      }                                 // Done with all partitions.
      // ------------------------------ //
      // Output the head M samples of the accumulator
      // ------------------------------ //
      for (size_t n=0;n<M;++n)          // For each sample in the block...
      {                                 // 
        outL[n]+=accL[n].real();        // Write the real part of the left channel accumulator to the output.
        outR[n]+=accR[n].real();        // Write the real part of the right channel accumulator to the output.
      }                                 // Done writing the output samples.
      // ------------------------------ //
      // Shift the accumulators by M samples and zero the tail
      // (OLA convolution).
      // ------------------------------ //
      std::rotate(accL.begin(),accL.begin()+M,accL.end()); // Shift the left channel accumulator by M samples.
      std::fill(accL.end()-M,accL.end(),std::complex<float>(0.0f,0.0f)); // Zero the tail of the left channel accumulator.
      std::rotate(accR.begin(),accR.begin()+M,accR.end()); // Shift the right channel accumulator by M samples.
      std::fill(accR.end()-M,accR.end(),std::complex<float>(0.0f,0.0f)); // Zero the tail of the right channel accumulator.
      // ------------------------------ //
      // Advance the partition index (if needed elsewhere).
      // ------------------------------ //
      if (++pidx>=part) pidx=0; // Advance the partition index and wrap it around if needed.
    }
    /// Convenience: process directly into AudioBuffer...
    void Process(const AudioBuffer<float>& in,AudioBuffer<float>& out) noexcept
    {
      assert(in.Channels()==1&&out.Channels()>=2);
      size_t N=in.Frames();
      auto* i0=in.Channel(0);
      auto* oL=out.Channel(0);
      auto* oR=out.Channel(1);
      ProcessBlock(i0,oL,oR);
    }
    /// API:
    double GetSampleRate(void) const noexcept { return fs; } // Get the sample rate.
    void SetSampleRate(double fs) noexcept
    {
        this->fs=fs;                     // Set the sample rate.
        engine.SetSampleRate(fs);        // Set the sample rate for the spectral operations engine.
    }
    size_t GetBlockSize(void) const noexcept { return M; } // Get the block size.
    size_t GetFFTSize(void) const noexcept { return fftsiz; } // Get the FFT size.
    size_t GetPartitionIndex(void) const noexcept { return pidx; } // Get the current partition index.
    float GetOverlap(void) const noexcept { return overlap; } // Get the overlap factor.
    Window<float>::WindowType GetWindowType(void) const noexcept
    {
      return win.GetWindowType();       // Get the current window type.
    }
    void SetWindowType(Window<float>::WindowType w) noexcept
    {
      win.SetWindowType(w,wSiz);        // Set the window type and size.
    }
    size_t GetIRLength(void) const noexcept { return irL.size(); } // Get the length of the impulse response.
    const auto& GetIRPartitions(void) const noexcept { return H_L; }
  private:
    size_t M=0;                         // Block size for the FFT convolution.
    size_t part=0;                      // Number of IR partitions.
    size_t fftsiz=32;                   // FFT size (must be power of 2).
    size_t pidx=0;                      // Current partition index.
    double fs=48000.0;                  // Sampling frequency (default 48kHz).
    std::vector<float> irL;             // Left channel impulse response.
    std::vector<float> irR;             // Right channel impulse response.
    Window<float> win;                  // Analysis window (2*M)
    SpectralOps<float> engine;          // Spectral operations engine for FFT/IFFT.
    float overlap=0.5f;                 // Overlap factor (default 50%).
    int wSiz=0;                         // Window size (2*M).
    std::vector<std::complex<float>> inBuf;// Input buffer for the FFT convolution.
    std::vector<std::complex<float>> accL;// OLA accumulator for the left channel.
    std::vector<std::complex<float>> accR;// OLA accumulator for the right channel.
    std::vector<std::vector<std::complex<float>>> H_L,H_R;// Partitioned FFTs
    size_t nextPowerOfTwo(size_t x) const noexcept
    {
      size_t n=1;                       // Start with first power of two.
      while (n<x) n<<=1;                // Shift left until we exceed x.
      return n;                         // Return the next power of two greater than or equal to x.;
    }                                   // -------- nextPowerOfTwo ---------------- //
    };
}

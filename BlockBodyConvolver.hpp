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
#include <cstring>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iterator>
#include "AudioBuffer.hpp"
#include "spectral/FCWTransforms.h"
#include "spectral/DSPWindows.h"
namespace dsp::wg
{
  using namespace dsp::spectral; // Use the spectral operations namespace.
  class BodyConvolverFFT
  {
  public:
    BodyConvolverFFT(void)=default; // Default constructor.
    void Prepare(double fs, size_t block) noexcept
    {                               // ---------- Prepare ------------------ //
      this->fs=fs;                // Set the sample rate.
      this->M=block;              // Set the block size.
      this->fftSize=nextPowerOfTwo(2*M); // Set the FFT size.
      // Scratch buffers for FFT/IFFT.
      xTime.assign(fftSize,{}); // Input block (time domain).
      xFreq.assign(fftSize,{}); // Input spectrum (frequency domain).
      yFreqL.assign(fftSize,{}); // Output spectrum for left channel.
      yFreqR.assign(fftSize,{}); // Output spectrum for right channel.
      // OLA output ring buffers.
      outOLAL.assign(M,0.0f);           // Overlap-add output ring buffer for left channel.
      outOLAR.assign(M,0.0f);           // Overlap-add output ring buffer for right channel.
      // Initialize the window object with Hanning window.
      win.GenerateWindow(Window<float>::WindowType::Hanning,2*M);
      // Shift by half-sample so w[0]=0.5,w[m]=0.5
      for (auto& w:win)
        w=std::sqrt(w);                 // sqrt(Hann)
      // A 50% overlap needs it centered  at M/2 but our window
      // is centered at M-1, so we rotate it.
      std::rotate(win.begin(),win.begin()+M/2,win.end());
      engine.SetSampleRate(fs);         // Set the sample rate for the spectral engine.
      engine.SetSamples(fftSize);       // Set the number of samples for the spectral engine.
      partCount=0;                      // Reset partition count.
      partIndex=0;                      // Reset partition index.
      prevIn.assign(M,0.0f);            // Previous input buffer for overlap-add.
    }                                   // ---------- Prepare ------------------ //
    bool LoadImpulseResponse(const std::vector<float>& L,const std::vector<float>& R)
    {
      if (L.empty()||L.size()!=R.size())
        return false;
      auto IRLen=L.size();              // Impulse response length.
      partCount=(IRLen+M-1)/M;          // The number of partitions needed.
      // Resize the partition buffers and circular buffer to fftSize length.
      H_L.assign(partCount,std::vector<std::complex<float>>(fftSize));
      H_R.assign(partCount,std::vector<std::complex<float>>(fftSize));
      Xq.assign(partCount,std::vector<std::complex<float>>(fftSize));
      // Local processing buffer.
      std::vector<std::complex<float>> buf(fftSize);
      // Process each partition of the impulse response.
      for (size_t p=0;p<partCount;++p)
      {
        // Clean previous impulse response residue
        // from the previous partition (what we are iterating over).
        std::fill(buf.begin(),buf.end(),std::complex<float>{});

        // Do the left channel first.
        for (size_t n=0;n<M;++n)
        {
          size_t idx=p*M+n;             // The idx in the partition.
          if (idx<IRLen)
            buf[n]={L[idx]*win[n],0.0f};// Apply the window's value at n to the impulse response.
        }                               // Done windowing the left channel.
        // Perform FFT on the left channel partition.
        H_L[p]=engine.FFT(buf);         // Compute the FFT of the left channel partition.
        // Clear the buffer for the right channel.
        std::fill(buf.begin(),buf.end(),std::complex<float>{});
        // And process the right channel.
        for (size_t n=0;n<M;++n)
        {
          size_t idx=p*M+n;             // The idx in the partition.
          if (idx<IRLen)
            buf[n]={R[idx]*win[n],0.0f};// Apply the window's value at n to the impulse response.
        }                               // Done windowing the right channel.
        // Perform FFT on the right channel partition.
        H_R[p]=engine.FFT(buf);         // Compute the FFT of the right channel partition.
      }                                 // Done processing all partitions.
      // Clear the circ buffer for input spectra.
      for (auto& v:Xq)
        std::fill(v.begin(),v.end(),std::complex<float>{});
      std::fill(outOLAL.begin(),outOLAL.end(),0.0f); // Clear the overlap-add output buffer for left channel.
      std::fill(outOLAR.begin(),outOLAR.end(),0.0f); // Clear the overlap-add output buffer for right channel.
      partIndex=0;                      // Reset the partition index.
      return true;                      // Return true to indicate success.
    }
  // ---------------------- Real time block ----------------------------------- //
  // Perfect-reconstruction WOLA (windowed overlap-add) convolution needs the same sqrt(Hann) at synthesis.
  // We do this by multiplying both havles (w0,w1) restores unity gain and tames
  // spectral leakage.
  void ProcessBlock(const float* in, float* outL, float* outR) noexcept
  {
    if (partCount==0)
    {
      std::memcpy(outL,in,M*sizeof(float)); 
      std::memcpy(outR,in,M*sizeof(float));
      return; // No IR loaded, just copy input to output.
    }
    if (partIndex==0)                       // Are we in the first FFT partition?
    {                                       // Yes, so we will...
      //std::fill(prevIn.begin(),prevIn.end(),in[0]);
      std::fill(prevIn.begin(),prevIn.end(),0.0f);// Start first window with silence.
    }                                      // Done checking if first FFT partition.
    assert (M&&"Prepare() not called"); // Ensure we have a block size.
    // 1. Window input into xTime and zero-pad.
    for (std::size_t n=0;n<M;++n)
    {                           // Feed 2*M samples to the FFT.
      // Previous half-window (tail of the last).
      xTime[n]={prevIn[n]*win[n],0.0f};       // Previous half-window samples.
      // Current half-window (head of the current).
      xTime[n+M]={in[n]*win[n+M],0.0f};        // Current half-window samples.
      // Store the current impulse sample because it'll be the 
      // prev of the next block.
      prevIn[n]=in[n];                        // Store the current input sample for the next block.
    }
    std::fprintf(stderr,"xTime0=%g\n",xTime[0].real()); // Debug: Print the first sample of xTime.
    // Zero-pad the rest (indices 2M ...ftSize-1).
    std::fill(xTime.begin()+2*M,xTime.end(),std::complex<float>{});
    // Perform FFT on xTime to get xFreq.
    //xFreq=engine.FFT(xTime);
    xFreq=engine.FFTStride(xTime);
    Xq[partIndex]=xFreq;          // Store the input spectrum in the circular buffer.
    // 2. Convolve with the IR spectra.
    std::fill(yFreqL.begin(),yFreqL.end(),std::complex<float>{}); // Clear the output spectrum for left channel.
    std::fill(yFreqR.begin(),yFreqR.end(),std::complex<float>{}); // Clear the output spectrum for right channel.
    for (size_t p=0;p<partCount;++p)
    {
      
      size_t q=(partIndex+partCount-p)%partCount; // Circular buf index.
      const auto& Xp=Xq[q];             // Get the input spectrum for this partition.
      const auto& HLp=H_L[p];           // Get the IR spectrum for left channel.
      const auto& HRp=H_R[p];           // Get the IR spectrum for right channel.
      // Convolve the input spectrum with the IR spectra.
      for (size_t k=0;k<fftSize;++k)
      {
        yFreqL[k]+=Xp[k]*HLp[k];        // Convolve left channel input with IR.
        yFreqR[k]+=Xp[k]*HRp[k];        // Convolve right channel input with IR.
      }                                 // Done convolving the input spectrum with the IR spectra.
    }                                   // Done processing all partitions.
    // 3. Perform IFFT on the output spectra.
    // Replaced FFT/IFFT with FFTStride/IFFTStride because it might be faster.
    // But FFTStride/IFFTStride already normalize so we comment out the for loops that normalized
    // after IFFT because we do not need them anymore.
    //xTime=engine.IFFT(yFreqL);          // Inverse FFT for left channel.
    //for (auto &c: xTime) c*=(1.f/fftSize);
    xTime=engine.IFFTStride(yFreqL);    // Inverse FFT the left channel.
    const auto yTimeL=xTime;            // Store the time domain output for left channel.
    //xTime=engine.IFFT(yFreqR);          // Inverse FFT for right channel.
    //for (auto &c:xTime) c*=(1.f/fftSize);
    xtime=engine.IFFTStride(yFreqR);    // Stride permutation IFFT on right channel.
    const auto yTimeR=xTime;            // Store the time domain output for right channel.
    // 4. Overlap-add the output.
    for (size_t n=0; n<M;++n)
    {
      const float w0=win[n];            // Sqrt(Hann) analysis window value at n.
      const float w1=win[n+M];          // ... second hald of the window (because win.size()==2).      float w
      // ------------------------------ //
      // Note: win[0…nZero-1] == 0 after the M/2 rotation, so the first few
      // samples of wetL / wetR are guaranteed to be 0.  This is expected –
      // they’ll be supplied by the overlap from the previous/next block.
      // ------------------------------ //
      float wetL=yTimeL[n].real()*w0+outOLAL[n];// Overlap-add for left channel.
      float wetR=yTimeR[n].real()*w0+outOLAR[n]; // Overlap-add for right channel.
      // Given that the OLA algorithm already folds the previous block's tail into outOLAL/outOLAR
      // We overwrite the value in outL/outR instead of accumulating it. Accumulating it would double
      // the history every window hop, which is not what we want.
      outL[n]=wetL;                     // Overwrite, not accumulate.
      outR[n]=wetR;                     // Overwrite, not accumulate.
      outOLAL[n]=yTimeL[n+M].real()*w1; // Stage next block overlap, windowed.
      outOLAR[n]=yTimeR[n+M].real()*w1; // Stage next block overlap, windowed.
    }                                   // Done OLA Processing.
    
    // 5. Update the partition index.
    partIndex=(partIndex+1)%partCount; // Increment the partition index and wrap around.
  }                                    // ------------ ProcessBlock --------------------- //
  void Process(const AudioBuffer<float>& monoIn,AudioBuffer<float>& stereoOut) noexcept
  {
    assert(monoIn.Channels()==1&&stereoOut.Channels()>=2);
    ProcessBlock(monoIn.Channel(0),stereoOut.Channel(0),stereoOut.Channel(1));
  }
  /// API:
  double GetSampleRate(void) const noexcept { return fs; } // Get the sample rate.
  size_t GetBlockSize(void) const noexcept { return M; } // Get the block size.
  size_t GetFFTSize(void) const noexcept { return fftSize; } // Get the FFT size.
  size_t GetPartitionIdx(void) const noexcept { return partIndex; } // Get the current partition index.
  private:
    // Global config:
    double fs=0.0;               // Sample rate.
    size_t M=0;                  // Block size.
    size_t fftSize=0;            // FFT size.
    size_t partCount=0;          // Number of partitions.
    size_t partIndex=0;          // Current partition index.
    Window<float> win;           // Window object for Hanning window.
    SpectralOps<float> engine;   // Spectral operations engine for FFT/IFFT.
    
    // IR Spectra and history.
    std::vector<std::vector<std::complex<float>>> H_L,H_R;

    // Input spectra circular buffer 
    std::vector<std::vector<std::complex<float>>> Xq;

    // Scratch buffers for FFT/IFFT.
    std::vector<std::complex<float>> xTime, xFreq;
    std::vector<std::complex<float>> yFreqL,yFreqR;

    // Overlap-add state
    std::vector<float> outOLAL, outOLAR; // Overlap-add output buffers for left and right channels.
    
    // Previous last M input samples
    std::vector<float> prevIn;
    size_t nextPowerOfTwo(size_t x) const noexcept
    {
      size_t n=1;                       // Start with first power of two.
      while (n<x) n<<=1;                // Shift left until we exceed x.
      return n;                         // Return the next power of two greater than or equal to x.;
    }                                   // -------- nextPowerOfTwo ---------------- //
    };
}

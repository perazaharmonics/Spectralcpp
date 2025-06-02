// ====================================================================================
// Description: A collection of common Spectral Transforms and manipulation methods.
//   used for various forms of Digital Signal Processing.
//
// Author: J. Enrique Peraza, JEP
// =======================================================================================

#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <vector>
#include <complex>
#include "DSPWindows.h"


using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Define alias templates for complex vectors and matrices.
/*template <typename T>
using complexV_t = vector<complex<T>>;
template <typename R>
using realV_t = vector<R>; // Template for real vectors
template <typename M>
using matrix_t = vector<vector<complex<M>>>; // Template for matrix vectors.*/

// Class template for spectral operations
template <typename T>
class Window;

template <typename T>
class SpectralOps
{
    using WindowType = typename Window<T>::WindowType; // Alias for WindowType
public:
    // Constructors
    SpectralOps(void);
    SpectralOps(const vector<T> &s);
    SpectralOps(const vector<T> &s, const double fs, const int len);
    SpectralOps(const vector<T> &s, const WindowType &w, const int windowSize);
    SpectralOps(const vector<T> &s, const WindowType &w, const int windowSize, const int overlap);
    ~SpectralOps(void);

    // Accessors
    const inline vector<T> GetSignal (void) const { return signal; }
    inline void SetSignal (const vector<complex<T>>& s) {singal=s;}
    const inline int GetSamples (void) const { return length; }
    inline void SetSamples (const int N) {length=N;}
    const inline double GetSampleRate (void) const { return sRate; }
    inline void SetSampleRate (const double fs) {sRate=fs;}
    const inline vector<complex<T>> GetTwiddles (void) const { return twiddles; }
    inline void SetSubCarrier const (const vector<complex<T>> &s) { subCarrier = s; }
    const inline vector<complex<T>> GetSubCarrier (void) { return subCarrier; } 
    
     // Operator Overloads
    vector<complex<T>> operator*(const vector<complex<T>> &s, const vector<complex<T>> &h);
    vector<complex<T>> operator*(const vector<complex<T>> &s, const vector<T> &h);
    vector<complex<T>> operator*(const vector<T> &s, const vector<complex<T>> &h);
    vector<complex<T>> operator*(const vector<T> &s, const vector<T> &h);
    vector<complex<T>> operator+(const vector<complex<T>> &s, const vector<complex<T>> &y);
   
    
    // Spectral Transformation methods.

    vector<complex<T>> FFTStride(const vector<complex<T>> &s);
    std::pair<vector<complex<T>>,vector<vector<complex<T>>>> FFTStrideEig(const vector<complex<T>> &s);
    vector<complex<T>> IFFTStride(const vector<complex<T>> &s);
    vector<complex<T>> FFT(const vector<complex<T>>& s);
    vector<complex<T>> IFFT(const vector<complex<T>>& s);
    vector<vector<complex<T>>> STFT(const vector<complex<T>> &s, const WindowType &w, const int wSiz, const float ovlap);
    vector<complex<T>> ISTFT(const vector<vector<complex<T>>> &sMat, const WindowType &w, const int wSiz, const float ovlap);
  // Spectral Computation Methods 
   
    vector<complex<T>> Convolution(const vector<complex<T>> &s, const vector<complex<T>> &h);        
    vector<complex<T>> OLAProcessor(const vector<complex<T>> &s,const WindowType &w, const int wSiz, const float ovlap);
    vector<complex<T>> OLAProcessor(const vector<complex<T>> &s, const vector<complex<T>> &h, const WindowType &w, const int wSiz, const float ovlap);
    vector<complex<T>> Shift(const vector<complex<T>>& s, const double fShift, const double fs);
    vector<vector<complex<T>>> Sweep (const vector<complex<T>>&s, const double fStart,const double fCenter, const double fStop,const double step,const double fs,const WindowType &w, const int wSiz, const float ovlap);
    // PSD should only be real to prevent aliasing and double.
    std::pair<std::vector<T>, T>LevinsonDurbin(const std::vector<T>& r, int order) const;
    std::vector<std::complex<T>> AR_PSD(const std::vector<T>& r, int order, int fftSize) const;
    vector<T> WelchPSD (const vector<T> &s, const WindowType& w, const int wSiz, const float ovlap, const int fftSiz);
protected:
  // Helper Methods
    int UpperLog2(const int N);
    void ForwardButterfly(vector<T> &last, vector<T> &curr, const vector<T> &twiddles, const int rot, const int nBits);
    void BitReversal(vector<T> &s, const int nBits);
    vector<T> ZeroPad(const vector<T> &s);
    vector<T> Upsample(const vector<T> &s, const int factor);
    vector<T> Downsample(const vector<T> &s, const int factor);
    vector<T> SetRBW(double rbw,double fs);
    vector<int> ToInt(const vector<complex<T>> &s);
    vector<double> ToReal (const vector<complex<T>> &s);    

private:
    vector<complex<T>> TwiddleFactor(int N); // Precompute the Twiddle for the FFT.
    vector<T> signal;                        // The signal to process.
    vector<T> subCarrier;                    // 
    WindowType window;                       // The window to apply to the signal.
    const int windowSize;                    // The size of the window.
    const float overlap;                     // The overlap factor.
    vector<complex<T>> twiddles;             // Precomputed twiddle factors.
    const double sRate;                      // The rate at which we sampled the RF.
    const int length;                        // The length of the signal.
};

// ------------------------------------ //
// Constructors and Destructors
// ------------------------------------ //
template <typename T>
SpectralOps<T>::SpectralOps(void) : signal{T(0)}, length{0}, sFreq{0.0}, window{WindowType::Rectangular}, windowSize{0}, overlap{0.5}
{
    // No window generation is needed for this constructor
}

template <typename T>
SpectralOps<T>::SpectralOps(const vector<T> &s) : signal{s}, length{0}, sFreq{0.0}, window{WindowType::Rectangular}, windowSize{0}, overlap{0.5}
{
    // No window generation is needed for this constructor
}

template <typename T>
SpectralOps<T>::SpectralOps(const vector<T> &s, const double fs, const int len) : signal{s}, length{len}, samplingFreq{fs}, window{WindowType::Rectangular}, windowSize{0}, overlap{0.5}
{
    // No window generation is needed for this constructor
}

template <typename T>
SpectralOps<T>::SpectralOps(const vector<T> &s, const WindowType &w, const int windowSize) : signal{s}, length{0}, samplingFreq{0.0}, window{w}, windowSize{windowSize}, overlap{0.5}
{
    // Generate the appropriate window for the given window size
    window = GenerateWindow(window, windowSize);
}

template <typename T>
SpectralOps<T>::SpectralOps(const vector<T> &s, const WindowType &w, const int windowSize, const int overlap) : signal{s}, length{0}, samplingFreq{0.0}, window{w}, windowSize{windowSize}, overlap{overlap}
{
    // Generate the appropriate window for the given window size
    window = GenerateWindow(window, windowSize);
}

template <typename T>
SpectralOps<T>::~SpectralOps(void)
{
    signal.clear();
    twiddles.clear();
}
// ======================= Operator Overload ================================ //
// Methods to manipulate the signals using C++ overloaded operators.
// ========================================================================== //
// Add two signals in the time domain.
template <typename T>
vector<complex<T>> operator+(const vector<complex<T>> &s, const vector<complex<T>> &y)
{
    vector<complex<T>> out;
    const T subCarrier = GetSubCarrier();      // Define the subcarrier frequency
    const T fs = GetSampleRate();      // Assume the sample rate corresponds to the size of the first signal

    if (s.size() != y.size()) 
    {
        cerr<<"Signals must be of the same size to add them."<<endl;
    }

    out.resize(s.size());

    for (size_t i = 0; i < s.size(); ++i) 
    {
      // Compute the complex exponential for the subcarrier modulation
      complex<T> modulator = polar(1.0, -2.0 * M_PI * subCarrier * i / fs);
        
      // Combine the two signals by adding the modulated subcarrier signal
      out[i] = s[i] + y[i] * modulator;
    }

    return out;
}


// ======================== Utility Methods ================================= //
// Utility Methods to precompute operations needed for Spectral Manipulations.
// ========================================================================== //
template <typename T>
vector<complex<T>> SpectralOps<T>::TwiddleFactor(int N)
{
    if (twiddles.size() != N / 2)                        // Did we precompute N/2 twiddle before?
    {                                                    // No, so we..
        twiddles.resize(N / 2);                          // Resize the twiddle factor vector.
        for (int i = 0; i < N / 2; ++i)                  //  loop for the N/2 points and
            twiddles[i] = polar(1.0, -2 * M_PI * i / N); //  compute the twiddle factors.
    }
    return twiddles;
}
// Get the smallest power of 2 that is greater than or equal to N
// that can hold the input sequence for the Cooley-Tukey FFT,
// which splits the input sequence into even and odd halves.
template <typename T>
int SpectralOps<T>::UpperLog2(const int N)
{
    for (int i = 0; i < 30; ++i) // For the first 30 powers of 2
      const int mask = 1 << i;   // Compute the value of 2^i
    if (mask >= N)               // If the power of 2 is >= N
        return i;                // Return the smallest power of 2 (i).
    return 30;                   // Else return 30 as the upper bound.
}

template <typename T>
vector<int>SpectralOps<T>::ToInt(const vector<complex<T>> &s)
{
    vector<int> sInt(s.size());
    for (size_t i = 0; i < s.size(); ++i)
        sInt[i] = static_cast<int>(s[i].real());
    return sInt;
}

template <typename T>
vector<double> SpectralOps<T>::ToReal(const vector<complex<T>> &s)
{
    vector<double> sReal(s.size());
    for (size_t i = 0; i < s.size(); ++i)
        sReal[i] = s[i].real();
    return sReal;
}
// Determine the amount of frequency bins to analyze per second of data.
template <typename T>
vector<T> SpectralOps<T>::SetRBW(double rbw, double fs)
{
  const int wSiz=static_cast<int>(fs/rbw);
  // Window is assumed to have been defined by the caller before calling this
  // method.
  return GenerateWindow(window,wSiz);
}

// ==================== Stride Permutation FFTs ============================= //
// Reference:  https://github.com/AndaOuyang/FFT/blob/main/fft.cpp
// ========================================================================== //
// Forward FFT Butterfly operation for the Cooley-Tukey FFT algorithm.
/* @param last: The previous stage of the FFT.
/*    Time domain signal iff first iteration of the FFT.
/*    Frequency domain signal iff IFFT.
/* @param curr: The temporary buffer for the FFT in this iteration.
/*   Frequency domain spectrum in the last iteration iff FFT
/*   Time domain signal in the last iteration iff IFFT.
/* @param twiddles: Vector of precomputed twiddle factors.
/* @param rot: The current stage of the FFT, iteration indicator. Starts at 0 for the first stage.
/* @param nBits: log2(N) where N is the length of the signal, total number of FFT stages.
/* Reference: https://github.com/AndaOuyang/FFT/blob/main/fft.cpp
*/
template <typename T>
void SpectralOps<T>::ForwardButterfly(vector<T> &last, vector<T> &curr, const vector<T> &twiddles, const int rot, const int nBits)
{
  if (rot == nBits)                          // Are we at the last stage of the FFT?
    return;                                  // Yes, so stop recursion.
    // ------------------------------------- //
    // Set the butterfuly section size to 2^(rot+1).
    // Each section doubles the size of the previous butterfly section.
    // ------------------------------------- //
  const int sectSiz = 1 << (rot + 1);        // Size of the butterfly section.
    // ------------------------------------- //
    // Number of sections (butterfly groups) the signal is split into at this stage. (phase groups) 
    // Each section is a group of butterflies, and has their phase computation.
    // ------------------------------------- //
  const int numSect = last.size() / sectSiz; // Number of sections the signal is divided into.
  const int phases = numSect;                // Number of phases (sections) in the FFT
    // ------------------------------------- //
    // Iterate over each phase in the FFT
    // Where each phase represents a group of butterfly operation
    // ------------------------------------- //
  for (int i = 0; i < phases; ++i)           // For every phase in the FFT
  {                                          // Perform the butterfly operation.
    const int base = i * sectSiz;            // Base index for the current phase.
    // ------------------------------------- //
    // Process each butterfly group within the current section.
    // The butterfly group is a pair of even and odd indices.
    // ------------------------------------- //
    for (int j = 0; j < sectSiz / 2; ++j) // For every butterfly group in the structure.
    {
    // ------------------------------------- //
    // Compute the even and odd indices in the butterfly group.
    // These elements will be combined to form the next stage of the FFT.
    // ------------------------------------- //      
        const int evenNdx = base + j;        // Even index in the butterfly group.
        const int oddNdx = base + sectSiz / 2 + j;// Odd index in the butterfly group.
    // ------------------------------------- //
    // Multiply the odd element by the twiddle factor for this butterfly group.  
    // The twiddle factor is a complex number that rotates the odd index.
    // and introduces the phase shift needed for the FFT. 
    // ------------------------------------- //   
        last[oddNdx] *= twiddles[j * phases];// Multiply the odd index by the twiddle factor.
    // ------------------------------------- //
    // Combine the next stage of the FFT using the even and odd indices.
    // The even and odd indices are combined to form the next stage of the FFT.
    // ------------------------------------- //      
        curr[evenNdx] = last[evenNdx] + last[oddNdx]; // Compute the even index.
        curr[oddNdx] = last[evenNdx] - last[oddNdx];  // Compute the odd index.
      } // Done with all butterfly groups.
  } // Done with all phases.
    // ------------------------------------- //
    // Recursivle move to the next stage of the FFT.
    // Swap the current and last buffers for the next iteration
    // ------------------------------------- //  
  ForwardButterfly(curr, last, twiddles, rot + 1, nBits); // Recurse to the next stage.
}
// Bit reversal permutation for the Cooley-Tukey FFT algorithm.
template <typename T>
void SpectralOps<T>:: BitReversal(vector<T> &s, const int nBits)
{
    // -------------------------------- //
    // Base Case: If the input size is <=2, no permutation necessary
    // For very small signals, bit reversal is not needed.
    // -------------------------------- //
  if (s.size()<=2)                      // Only two or less samples?
    return;                             // Yes, so no need to reverse bits.
    // -------------------------------- //
    // Special Case: If the input is exactly 4 samples, swap the middle
    // two elements. Handle the 2-bit case directly.
    // -------------------------------- //
  if (s.size()==4)                      // Is the signal exactly 4 samples?
  {                                     // Yes, so swap the middle two elements.
    swap(s[1], s[2]);                   // Swap the middle two elements.
    return;                             // Done with the bit reversal.
  }
    // -------------------------------- //
    // General Case: For signals larger than 4 samples, perform bit reversal.
    // Initialize a vector to hold bit-reversed indices and compute the bit
    // reversed indices for the FFT.
    // -------------------------------- //
  vector<int> revNdx(s.size());         // Vector to hold bit-reversed indices.
    // -------------------------------- //
    // Manually set the first 4 indices' bit-reversed values.
    // These are the known bit reversed values for the 2-bit case.
    // -------------------------------- //
  revNdx[0]=0;                          // Bit-reversed index for 0 is 0.
  revNdx[1]=1<<(nBits-1);               // == 100...0 in binary == 2^(nBits-1).
  revNdx[2]=1<<(nBits-2);               // == 010...0 in binary == 2^(nBits-2).
  revNdx[3]=revNdx[1]+revNdx[2];        // == 110...0 in binary == 2^(nBits-1) + 2^(nBits-2).
    // -------------------------------- //
    // Loop through to  compute the rest of the bit-reversed indices.
    // the bit-reversed index is the reverse of the binary representation of the index.
    // -------------------------------- //
    // Theorem: For all nk=2^k-1 where k<= nBits, 
    // revNdx[nk]=revNdx[n(k-1)]+2^(nBits-k)
    // revNdx[nk-i]=revNdx[nk]-revNdx[i]
    // -------------------------------- //
  for (int k=3; k<=nBits;++k)           // For all remaining bits in the signal.
  {
    const int nk=(1<<k)-1;              // Compute nk=2^k-1.
    const int nkmin1=(1<<(k-1))-1;      // Compute n(k-1)=2^(k-1)-1.
    // -------------------------------- //
    // Derive the bit-reversed index for nk using the bit reversal of n(k-1).
    // The bit-reversed index for nk is the bit-reversed index for n(k-1) plus 2^(nBits-k).
    // -------------------------------- //
    revNdx[nk]=revNdx[nkmin1]+(1<<(nBits-k)); // Compute revNdx[nk].
    // -------------------------------- //
    // Loop to compute the remaining bit reversed indices.
    // Compute for the range nk -i using nk and previously computed values.
    // -------------------------------- //
    for (int i=1; i<=nkmin1;++i)        // For the range nk-i.
      revNdx[nk-i]=revNdx[nk]-revNdx[i]; // Compute revNdx[nk-i].
  }
    // -------------------------------- //
    // Permute the signal using the bit-reversed indices.
    // Swap elements if the bit-reversed index is greater than the current index.
    //--------------------------------- //
  for (int i=0; i<s.size();++i)         // For all elements in the signal.
    if (i<revNdx[i])                    // If the index is less than the bit-reversed index.
      swap(s[i], s[revNdx[i]]);         // Swap the elements.             
}                                       // End of the function.
    // ------------------------------------------------------------------------
    // LevinsonDurbin: Given autocorrelation r[0..p], solves for AR coefficients
    //   r[0] a[0] + r[1] a[1] + … + r[p] a[p] = 0,    (Toeplitz system)
    //   returns (a[1..p], σ²) where σ² is the final prediction error.
    //   “order” = p.  We assume r.size() >= p+1.
    // ------------------------------------------------------------------------
    template <typename T>
    std::pair<std::vector<T>, T>SpectralOps<T>::
    LevinsonDurbin(const std::vector<T>& r, int order) const
    {
        // r: autocorrelation, r[0] … r[order]
        // order: AR order (p)
        if ((int)r.size() < order+1) {
            throw std::invalid_argument{"LevinsonDurbin: need r.size() >= order+1"};
        }
        std::vector<T> a(order+1, T{0}); // a[0]..a[p], we keep a[0]=1 internally
        std::vector<T> e(order+1, T{0}); // prediction error at each stage
        a[0] = T{1};
        e[0] = r[0];
        if (std::abs(e[0]) < std::numeric_limits<T>::epsilon()) {
            // All‐zero autocorrelation → trivial
            return { std::vector<T>(order, T{0}), T{0} };
        }

        for (int m = 1; m <= order; ++m) {
            // Compute reflection coefficient κ_m
            T num = r[m];                  // numerator = r[m] + sum_{i=1..m−1} a[i]·r[m−i]
            for (int i = 1; i < m; ++i) {
                num += a[i] * r[m - i];
            }
            T kappa = - num / e[m-1];

            // Update a[1..m]:
            std::vector<T> a_prev(m+1);
            for (int i = 0; i <= m; ++i) a_prev[i] = a[i];
            a[m] = kappa;
            for (int i = 1; i < m; ++i) {
                a[i] = a_prev[i] + kappa * a_prev[m - i];
            }

            // Update prediction error
            e[m] = e[m-1] * ( T{1} - kappa * kappa );
            if (std::abs(e[m]) < T{0}) {
                e[m] = T{0};
            }
        }

        // Return only a[1..p] (drop a[0]=1) and final error e[p]
        std::vector<T> arCoeffs(order);
        for (int i = 1; i <= order; ++i) {
            arCoeffs[i-1] = a[i];
        }
        return { arCoeffs, e[order] };
    }

    // ------------------------------------------------------------------------
    // AR_PSD: Given autocorrelation r[0..p], compute the “all‐pole” PSD estimate
    //    at fftSize uniformly spaced frequencies [0, 2π).  We solve AR(p) via
    //    Levinson‐Durbin, then evaluate
    //      H(ω) = σ² / |1 + a[1] e^{-jω} + … + a[p] e^{-j p ω} |²
    //    at Nfft points, returning a vector<complex<T>> of length Nfft
    //    (you can take real(H) or abs(H)² as your PSD). 
    // ------------------------------------------------------------------------
    template <typename T>
    std::vector<std::complex<T>>SpectralsOps<T>::
    AR_PSD(const std::vector<T>& r, int order, int fftSize) const
    {
        if (order < 1 || (int)r.size() < order+1) {
            throw std::invalid_argument{"AR_PSD: order must be ≥1 and r.size() ≥ order+1"};
        }
        // 1) run Levinson‐Durbin on r[0..order]
        auto [a, sigma2] = LevinsonDurbin(r, order);
        // a = vector length p, contains a[1],…a[p], and sigma2 = error at final stage

        // 2) build PSD at fftSize freq bins
        std::vector<std::complex<T>> psd(fftSize);
        const T normFactor = T{2} * M_PI / static_cast<T>(fftSize);
        for (int k = 0; k < fftSize; ++k) {
            T omega = normFactor * static_cast<T>(k); 
            // Evaluate denominator D(ω) = 1 + ∑_{m=1..p} a[m-1] e^{-j m ω}
            std::complex<T> denom = T{1};
            for (int m = 1; m <= order; ++m) {
                denom += a[m-1] * std::exp(std::complex<T>(T{0}, -omega * static_cast<T>(m)));
            }
            // PSD(ω_k) = σ² / |D(ω)|²
            std::complex<T> H = std::complex<T>(sigma2) / (denom * std::conj(denom));
            psd[k] = H;
        }
        return psd;
    }

// ======================== Stride FFTs ===================================== //
// Stride FFTs are a special case of the FFT that uses a stride to compute the FFT
// of a signal. This is useful for signals that are not a power of 2 in length.
// The FFTStride method computes the FFT of a signal using the Cooley-Tukey algorithm
// with a stride. The IFFTStride method computes the IFFT of a signal using the
// FFTStride method with the conjugate of the input signal.
// ========================================================================== //
// FFTStrideEig computes the FFT of a signal and returns the spectrum and the
// eigenvectors of the FFT matrix. This is useful for computing the eigenvalues
// and eigenvectors of the FFT matrix, which can be used for spectral analysis
// to obtain phase information.
template <typename T>
std::pair<vector<complex<T>>,vector<vector<complex<T>>>> SpectralOps<T>::FFTStrideEig(const vector<complex<T>> &s)
{
  if (s.empty())                        // Is the input signal empty?
    return {vector<complex<T>>(), vector<vector<complex<T>>>()}; // Yes, so return empty vectors.
  // ---------------------------------- //
  // Calculate the number of bits needed for the FFT rounded to the
  // nearest upper power of 2. This is the number of stages in the FFT butterfly.
  // ---------------------------------- //
  const int NBits=UpperLog2(static_cast<int>(s.size())); // Get the number of bits for the FFT.
  const int N=1<<NBits;                 // Get the FFT length as a power of 2.
  // ---------------------------------- //
  // Create temporary buffers for the FFT.
  // The last buffer holds the previous stage of the FFT.
  // The current buffer holds the current stage of the FFT.
  // ---------------------------------- //
  vector<complex<T>> last(N), curr(N);  // Temporary buffers for the FFT.
  // ---------------------------------- //
  // Copy the input signal to the last buffer, and zero-pad if necessary.
  // ---------------------------------- //
  copy(s.begin(), s.end(), last.begin()); // Copy the input signal to the last buffer.
  // Zero-pad the last buffer to the FFT length.
  if (s.size() < N)                     // Is the input signal smaller than the FFT length?
    fill(last.begin() + s.size(), last.end(), complex<T>(0)); // Yes, so zero-pad the last buffer.
  // ---------------------------------- //
  // Perform the bit reversal permutation on the input signal.
  // This reorders the input signal to prepare for the Cooley-Tukey FFT.
  // ---------------------------------- //
  BitReversePermutation(last, NBits);   // Perform bit reversal permutation.
  // ---------------------------------- //
  // Perform the FFT butterfly operation for the Cooley-Tukey FFT.
  // This computes the FFT in-place using the last and current buffers.
  // This is where the Cooley-Tukey FFT algorithm takes place.
  // ---------------------------------- //
  ForwardButterfly(last, curr, twiddles, 0, NBits); // Perform the FFT butterfly.
  // ---------------------------------- //
  // Here we compute the Fourier matrix and index the eigenvectors.
  // Return the computed FFT spectrum and the eigenvectors of the FFT matrix.
  // The eigenvectors are the Fourier basis vectors, which are the twiddle factors.
  // ---------------------------------- //
  vector<vector<complex<T>>> eigvecs(N, vector<complex<T>>(N)); // Create a matrix for the eigenvectors.
  const T invsqrt=static_cast<T>(1)/sqrt(static_cast<T>(N)); // Inverse square root of N for normalization.
  for (int k=0;k<N;k++)                // For each row in the eigenvector matrix...
  {                                    // Compute the Fourier matrix.
    long double angle=-2.0L*M_PI*std::static_cast<long double>(ell)*std::static_cast<long double>(k)/std::static_cast<long double>(N);
    eigvecs[ell][k]=complex<T>(std::cos(angle),std::sin(angle))*invsqrt; // Compute the k-th eigenvector.
  }
}

template <typename T>
vector<complex<T>> SpectralOps<T>::FFTStride (const vector<complex<T>> &s)
{
    // ---------------------------------- //
    // Base Case: If the input is empty, return an empty vector.
    // ---------------------------------- //
    if (s.empty())                        // Is the input signal empty?
        return vector<complex<T>>();      // Yes, so return an empty vector.
    // ---------------------------------- //
    // Calculate the number of bits needed for the FFT rounded to the 
    // nearest upper power of 2. This is the number of stages in the FFT.
    // ---------------------------------- //
    const int nBits=UpperLog2(s.size());  // Get the number of bits for the FFT.
    // ---------------------------------- //
    // Calculate the FFT length as 2^nBits.
    // This is the length of the FFT signal.
    // ---------------------------------- //
    const int N=1<<nBits;                 // Get the FFT length as a power of 2.
    // ---------------------------------- //
    // Precompute the twiddle factors for the FFT.
    // The twiddle factors are used to rotate the signal in the FFT.
    // ---------------------------------- //
    const vector<complex<T>> twiddles=TwiddleFactor(N); // Phase-frequency vector.
    // ---------------------------------- //
    // Create temporary buffers for the FFT.
    // The last buffer holds the previous stage of the FFT.
    // The current buffer holds the current stage of the FFT.
    // ---------------------------------- //
    vector<complex<T>> last(N), curr(N);  // Temporary buffers for the FFT.
    // ---------------------------------- //
    // Copy the input signal to the last buffer, and zero-pad if necessary.
    // ---------------------------------- //
    copy(s.begin(), s.end(), last.begin()); // Copy the input signal to the last buffer.
    // ---------------------------------- //
    // Perform the bit reversal permutation on the input signal.
    // This reorders the input signal to prepare for the Cooley-Tukey FFT.
    // ---------------------------------- //
    BitReversePermutation(last, nBits);   // Perform bit reversal permutation.
    // ---------------------------------- //
    // Perform the FFT butterfly operation for the Cooley-Tukey FFT.
    // This computes the FFT in-place using the last and current buffers.
    // This is where the Cooley-Tukey FFT algorithm takes place.
    // ---------------------------------- //
    ForwardButterfly(last, curr, twiddles, 0, nBits); // Perform the FFT butterfly.
    // ---------------------------------- //
    // Return the computed FFT spectrum.
    // ---------------------------------- //
    if (nBits %2 == 1)                    // Is the number of bits odd?
        return curr;                      // Yes, so return the current buffer.
    return last;                          // No, so return the last buffer.
}
template <typename T>
// The IFFT can be computed using the FFT with flipped order of the 
// frequency bins. That is, the complex conjugate of the input signal.
//   and thus the twiddle factors.
// So we just flip the frequency spectrum an normalize by 1/N.
// ------------------------------------------------------------
// Theorem: Let x[n] denote a time-domain signal and X[k] denote a frequency
// domain signal,then: 
// x[n]=(1/N) * SUM[k=0 to N-1] * {X[k] * exp(j*(2*pi/N)*k*n)} == IFFT(X[k]) 
// Let's denote m=-k, then: 
// x[n]=(1/N)*SUM[m=0 to 1-N]*{X[m]*exp(-j*(2*pi/N)*k*n)==FFT(X[m])
// We know that FFT is circularly periodic, thus X[m]=X[-k]=X[n-k]. 
// Therefore we can get X[m], simply by reversing the order of X[k].
// --------------------------------------------------------------
vector<complex<T>> SpectralOps<T>:: IFFTStride (const vector<complex<T>>& s)
{
  vector<complex<T>> sConj(s);          // Copy the input signal.
  // ---------------------------------- //
  // Flip the frequency spectrum
  // ---------------------------------- //
  reverse(next(sConj.begin()),sConj.end()); // Reverse the frequency spectrum.
  const double siz=sConj.size();        // The length of conjugated spectrum.
  // ---------------------------------- //
  // Normalize the signal by 1/N using lambda function.
  // ---------------------------------- //
  transform(sConj.begin(), sConj.end(), sConj.begin(), 
    [siz](complex<T> x){return x/static_cast<T>(siz);}); // Normalize the signal.
  return FFTStride(sConj);              // Return the FFT of the conjugate.
}

// ============================= FFT and IFFT Algorithms ============================= //

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// So to get the FFT of a signal x(n) of length N we have to divide and conquer.
// We do this by using the Cooley-Tukey Algorithm, described here as:
// 1. Divide the signal into even and odd samples, so we have:
//      for (i.begin(); i.end() )
//        evenT[i]=x[i*2]
//        oddT[i]=x[i*2+1]
// 2. Conquer the signal; we recursively apply the FFT on both halves:
//      X_even(k) = FFT(evenT)
//      X_odd(k) = FFT(oddT)
// 3. Now we precompute the twiddle factors, this saves A LOT of time:
//      TwiddleFactor(k) = exp( -j * (2*pi/n)*k)
// 4. Having the Twiddle Factors we compute the FFT butterfly to obtain the full
//    frequency spectrum - the amplitude and phase of sines and cosines that
//    composit it.
//      for (k.begin; k.end/2)
//        t=TwiddleFactor(k)*X_odd(k)
//        X(k)=X_even(k)+t
//        X(k+N/2)=X_even(k)-t
// 5. Return the spectrum of the signal
// Note that this algorithm should be slower than the FFT Stride above, but it 
// is also clearer. 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename T>
vector<complex<T>> SpectralOps<T>:: FFT(const vector<complex<T>>& s)
{
  const int N=s.size();                 // The length of the input signal.
    // -------------------------------- //
    // Base Case: When the input is 1, return the signal.
    // -------------------------------- //
  if (N<=1)                             // Is it a single point?
    return s;                           // Yes, return the signal.
    // -------------------------------- //
    // Divide Step: Divide the signal into even and odd samples.
    // -------------------------------- //
  vector<complex<T>> evenT(N/2), oddT(N/2);
  for (int i=0; i<N/2;++i)              // Up to the folding frequency.
  {
    evenT[i]=s[i*2];                    // Get the even samples.
    oddT[i]=s[i*2+1];                   // Get the odd samples.
  }                                     // Done decimating in time.
    // -------------------------------- //
    // Conquer Step: Recurse, apply FFT to evens and odds.
    // -------------------------------- //
  vector<complex<T>> evenF=FFT(evenT);  // Transform even samples.
  vector<complex<T>> oddF=FFT(oddT);    // Transform odd samples.
    // -------------------------------- //
    // Precompute the twiddle factors.
    // -------------------------------- //
  vector<complex<T>> tf=TwiddleFactor(N);// Get the phase-freq rotation vector.
    // -------------------------------- //
    // Compute the FFT butterfly for this section
    // -------------------------------- //
  vector<complex<T>> S(N);              // Initialize freq-domain vector.
  complex<T> t{0.0,0.0};                // Single root of unity.
  for (int k=0; k<N/2; ++k)             // Up to the folding frequency.
  {
    // -------------------------------- //
    // Get the amplitude phase contribution for current butterfly phase.
    // -------------------------------- //
    t=twiddle[k]*oddF[k];               // Scale and get this freq bin.
    // -------------------------------- //
    // Prepare results for next butterfly phase.
    // -------------------------------- //
    S[k]=evenF[k]+t;                    // Produce even result for nxt butterfly.
    S[k]=oddF[k]-t;                     // Produce odd result for nxt butterfly.
  }                                     // Done computing butterfly.
    // -------------------------------- //
    // Return computed spectrum. 
    // -------------------------------- //
  return S;                             // The computed spectrum.
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ IFFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// So how do we calculate the Inverse FFT? To get the inverse of a signal X(k)
// of length k an elegant trick is performed.
// The IFFT is nothing more than the FFT multiplied by 1/N, and with a twiddle
// factor that rotates clockwise, instead of counter-clockwise.
// It is nothing more than the conjugate of the FFT multiplied by (1/N).
//
// The operation goes as follows:
// 1. Conjugate the discrete frequency input signal X(k):
//    X_conjugate(k) = Re(X(k)) - j*Im(X(k))
// 2. Next we perform the FFT on the conjugated signal, this performs the IDFT
//    but does not scale - that comes afterwards:
//    X(k) = SUM[n=0 to N-1] {x(n) * exp(-j*(2*pi/N)*k*n) }
// 3. Now we conjugate again and multiply by (1/N), this returns our signal to
//    the time domain: (wizardry at hand!)
//    x(n) = (1/N) * SUM[k=0 to N-1] * {X(k) * exp(-j*(2*pi/N)*k*n)}
// 4. Return the time-domain samples of the signal.
// Note that this algorithm should be slower than the IFFT Stride above, but it 
// is also clearer.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename T>
vector<complex<T>> SpectralOps<T>:: IFFT (const vector<complex<T>> &s)
{
  const int N=s.size();                 // Get the length of the signal.
    // -------------------------------- //
    // 1. Reverse the frequency spectrum by conjugating the input signal.
    // -------------------------------- //
  vector<complex<T>> sConj(N);          // Reversed spectrum buffer.
  for (int i=0; i<N;++i)                // For every sample in the spectrum
    sConj[i]=conj(s[i]);                //   reverse the frequency bins.
    // -------------------------------- //
    // 2. Perform FFT on the conjugated spectrum.
    // -------------------------------- //
  vector<complex<T>> S=FFT(sConj);      // Reverse-spectrum buffer.
    // -------------------------------- //
    // 3. Conjugate and normalize the signal.
    // -------------------------------- //
  vector<complex<T>> s(N);              // Signal buffer.
  for (int i=0;i<N;++i)                 // For all samples in reversed spectrum
    s[i]=conj(S[i])/static_cast<T>(N);  // Reverse and normalize.
  return s;                             // Return time-domain signal.  
} 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Convolution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Method to perform the convolution of two signals. Typically in our context
// the convolution is done between an input signal s(n) of length N  and
// filter of length N h(n).
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename T>
vector<complex<T>> SpectralOps<T>::Convolution( // Obtain the product of two signals.
    const vector<complex<T>> &s,        // Our input signal.
    const vector<complex<T>> &h)        // Our filter.
{                                       // ---------- Convolution ----------- //
    const int n = s.size();             // Length of the input signal.
    const nt m = h.size();              // Length of the filter.
    int N = 1;                          // Size of the N-point FFT.
    while (N < n + m - 1)               // While N is less than the convolution
        N <<= 1;                        // Find smallest power of 2 >= n+m-1
    // -------------------------------- //
    // Zero-Pad the signal and filter to the length of the N-Point FFT.
    // -------------------------------- //
    vector<complex<T>> sPad = ZeroPad(s); // Zero pad the input signal.
    vector<complex<T>> hPad = ZeroPad(h); // Zero pad the filter.
    // -------------------------------- //
    // Apply the FFT on the Zero-Padded signal.
    // -------------------------------- //
    vector<complex<T>> S = FFT(sPad);   // The FFT of the input signal.
    vector<complex<T>> H = FFT(hPad);   // The FFT of the filter.
    // -------------------------------- //
    // Now the filtered signal is just the product of their spectrum.
    // -------------------------------- //
    vector<complex<T>> Y(N);            // Place to store resulting spectrum.
    for (int i = 0; i < N; ++i)         // For the N-Points of the FFT.
        Y[i] = S[i] * H[i];             // Get the filtered spectrum.
    // -------------------------------- //
    // Obtain the time-domain resulting signal using the IFFT.
    // -------------------------------- //
    vector<complex<T>> y(N);            // Place to store resulting signal.
    y = IFFT(Y);                        // Get the resulting signal.
    y.resize(n + m - 1);                // Truncate to the original size.
    return y;                           // Return the filtered signal.
}                                       // ---------- Convolution ----------- //
// Short-Time Fourier Transform
template <typename T>
vector<vector<complex<T>>> SpectralOps<T>::STFT(const vector<complex<T>> &s, 
  const WindowType &w, 
  int wSiz, 
  const float overlap)
{
  int step = wSiz * (1 - overlap / 100.0);// step size based on the ovlap %
  int nSegs = (s.size() - wSiz + step) / step;// # segs for the signal.
  vector<vector<complex<T>>> sMat(nSegs, vector<complex<T>>(wSiz));
  vector<T> window = GenerateWindow(windowType, wSiz);// Window to be applied.
    // -------------------------------- //
    // Process each seg of the signal. Each row of the matrix is a frame of
    // the signal, and every column a frequency bin inside the windowed segment.
    // -------------------------------- //
  for (int i=0;i<nSegs;++i)
  {   
      int start = i * step;             // Starting ndx for our current segment.
      vector<complex<T>> seg(wSiz, 0.0); // Initialize the seg with zeros
    // -------------------------------- //    
    // Apply the window to the segment and copy the windowed signal.
    // For the size of the window and remaining frequency bins of the signal.
    // -------------------------------- //
      for (int j = 0; j <wSiz && start+j< s.size(); ++j)
        seg[j] = s[start + j] * window[j];
    // -------------------------------- //
    // Compute the FFT of the windowed seg and store it in the STFT matrix
    // -------------------------------- //
      sMat[i] = FFT(seg, windowType);
  }
  return sMat;                          // The windowed spectrum.
}
template <typename T>
// Inverse Short-Time-Fourier-Transform Method.
vector<complex<T>> SpectralOps<T>::ISTFT(
    const vector<vector<complex<T>>> &sMat,// Input STFT matrix
    const WindowType &w,                // Window type
    const int wSiz,                     // Window size (should be power of 2)
    const float ovlap)                  // Overlap percentage (e.g., 50%)
{
    int step = wSiz*(1 -ovlap/100.0);// step size based on the ovlap %
    int len = step*(sMat.size()-1)+wSiz;// Length of the original signal.
    // --------------------------------- //
    // Calculate the overlap count for the normalization step.
    // -------------------------------- //
    vector<complex<T>> sig(len,0.0);// Initialize the result signal
    vector<T> nOverlaps(len,0.0); // Initialize the overlap count
    // Generate the window to be applied during the inverse process
    vector<T> window = GenerateWindow(w, wSiz);
    // -------------------------------- //
    //  Process each seg of the signal. Each row of the matrix is a frame of
    //  of the signal, and each column a frequency bin inside the window.
    // -------------------------------- //
    for (int i = 0; i < sMat.size(); ++i)
    {
      int start = i*step;               // Starting ndx of the frame (segment).
    // -------------------------------- //
    // Compute the IFFT of the current segment, and get the time-domain short
    // signal. This allows us to reconstruct it back using its segments. 
    // -------------------------------- //
      vector<complex<T>> seg = IFFT(sMat[i]);
    // -------------------------------- //
    // Overlap-add the IFFT result to the output signal. Because the segments
    // were windowed, we overlap-add each segment to obtain total signal's energy
    // contribution. Because the energy of the signal is greater in the
    // samples that lie in the overlapped region, we keep track of these in the 
    // second step. This allows us not to overshoot the signal's original ampli
    // -tude in these regions, as the signal's energy is purposely overlapped
    // by the spectral windows.
    // -------------------------------- // 
      for (int j=0;j<wSiz && (start+j)<sig.size();++j)
      {
        sig[start+j]+=seg[j]*window[j]; // Frame contribution to long signal.
        nOverlaps[start+j]+=window[j];  // Keep track of Overlapp-Add Windowing process 
      }
    }
    // -------------------------------- //
    // Normalize the result by dividing by the overlap count. We want to know 
    // how many times each segment of the signal was affected by the OLA Window
    // process, so that we can scale the original signal's amplitude accurately.
    // -------------------------------- //
    for (int i = 0; i < sig.size(); ++i)
      if (nOverlaps[i] != 0.0)
        sig[i] /= nOverlaps[i];
    // Return the reconstructed time-domain signal
    return sig;
}
// vector<complex<T>> OLAProcessor(const vector<complex<T>> &s, const vector<complex<T>> &h, const WindowType &w, const int wSiz, const float ovlap)
template <typename T>
vector<complex<T>> SpectralOps<T>::OLAProcessor(
    const vector<complex<T>> &s, // The input signal.
    const vector<complex<T>> &h, // The desired FIR filter.
    const WindowType &w,         // The window used.
    const int wSiz,              // The size of the window.
    const float ovlap)           // The percentage of overlap
{
    vector<vector<complex<T>>> sMat = STFT(s, w, wSiz, ovlap); // STFT of input signal.
    vector<complex<T>> H = FFT(h);      // FFT of the FIR filter.
    const int frames = sMat.size();     // Number of frames in the STFT matrix.
    vector<vector<complex<T>>> sig(frames, vector<complex<T>>(wSiz));
    // -------------------------------- //
    // Perform element-wise multiplication of the STFTs
    // -------------------------------- //
    for (int i = 0;i<frames;++i)        // For the number of frames...
      for (int j=0;j <wSiz; ++j)        // .. and the freq bins in the window.
        sig[i][j]=sMat[i][j]*H[j];      // Convlute the short signal.
    // -------------------------------- //
    // Perform the inverse STFT to get the filtered time-domain signal.
    // -------------------------------- //
    return ISTFT(sig, w, wSiz, ovlap);
}
// vector<complex<T>> OLAProcessor(const vector<complex<T>> &s, const vector<complex<T>> &h, const WindowType &w, const int wSiz, const float ovlap)
template <typename T>
vector<complex<T>> SpectralOps<T>::OLAProcessor(
    const vector<complex<T>> &s, // The input signal.
    const WindowType &w,         // The window used.
    const int wSiz,              // The size of the window.
    const float ovlap)           // The percentage of overlap
{
    vector<vector<complex<T>>> sMat = STFT(s, w, wSiz, ovlap);// STFT of input signal.
    vector<complex<T>> H(wSiz,complex<T>(0.0,0.0));// Dummy Impulse filter.
    vector<complex<T>> H = FFT(h);      // FFT of the FIR filter.
    const int frames = sMat.size();     // Number of frames in the STFT matrix.
    vector<vector<complex<T>>> sig(frames, vector<complex<T>>(wSiz));
    // -------------------------------- //
    // Perform element-wise multiplication of the STFTs
    // -------------------------------- //
    for (int i = 0;i<frames;++i)        // For the number of frames...
      for (int j=0;j <wSiz; ++j)        // .. and the freq bins in the window.
        sig[i][j]=sMat[i][j]*H[j];      // Convlute the short signal.
    // -------------------------------- //
    // Perform the inverse STFT to get the filtered time-domain signal.
    // -------------------------------- //
    return ISTFT(sig, w, wSiz, ovlap);
}
// Determine the power sepctral density of a windowed signal using Welch's method.
template <typename T>
vector<T> SpectralOps<T>::WelchPSD(
 const vector<T> &s,                    // The signal to process.
 const WindowType& w,                   // The window to apply to the signal.
 const int wSiz,                        // The size of the window.
 const float ovlap,                     // Overlap percentage (50% typical)
 const int fftSiz)                      // The size of the FFT.
{
    // -------------------------------- //
    // Compute the STFT of the signal.
    // -------------------------------- //
  vector<vector<complex<T>>> stftMat=STFT(signal,wType,wSiz,ovlap);
    // -------------------------------- //
    // Determine the scale factor of the window.
    // -------------------------------- //
  vector<T> pxxAvg(N,T(0));             // Initialize PSD Buffer
  const double winScl=pow(norm(GenerateWindow(w,wSiz),2),2); // Get Scale Factor.
    // -------------------------------- //
    // Now we accumulate the PSD for each segment in the STFT matrix.
    // -------------------------------- //
  for (int i=0; i<stftMat.size();++i)   // For each fft segment.
  {                                     // Where each row in the STFT matrix is
    const vector<complex<T>>& fftSegment=stftMat[i];//  an FFT segment
    // -------------------------------- //
    // Compute the Power Spectal Density
    // -------------------------------- //
    for (int j=0;j<fftSiz;++j)          // For every sample to process by the FFT
      pxxAvg[j]+=norm(fftSegment[j])/winScl;// Accumulate sq. magnitude.
  }                                     // Done accumulating sq.magntiude segments.
    // -------------------------------- //
    // Now we average the PSD over the segments
    // -------------------------------- //
  for (int i=0; i<fftsiz; ++i)          // For every sample to process by the FFT
    pxxAvg[i]/=stftMat.size();          // Average PSD over all segments.
    // -------------------------------- //
    // and normalize the PSD over 2*pi periodic interval.
    // -------------------------------- //
  for (int i=0; i<fftSiz;++i)           // For every sample...
    pxxAvg[i]/=(2*M_PI);                //
    // -------------------------------- //
    // Ensure the total energy of the sigal is conserved in the sprectrm.
    // Parseval's Theorem: SUM[n to N-1] x[n]^2 = (1/N) SUM[n to N-1] X[k]^2
    // -------------------------------- //
  pxxAvg[0]/=2;                         // Avergave freq bin 1 (DC component).
  for (int i=0;i<fftSiz;++i)            // For ever sample count the energy in
    pxxAvg[i]*=2;                       //  positive and negative halves
    // -------------------------------- //
    // Return the first half of the PSD (our FFT is symmetric) and we already
    // have recollected all the power.
    // -------------------------------- //
  return vector<T>(pxxAvg.begin(),pxxAvg.being()+fftSiz/2+1);

}
template <typename T>
// Method to perform a frequency shift of the center frequency by a const amount
vector<complex<T>> SpectralOps<T>::Shift(  // Shift the signal in frequency domain.
  const vector<complex<T>>& s,           // The input signal.
  const double fShift,                  // The amount to shift it by
  const double fs)                      // The sample rate of the signal
{                                       // ---------- Shift ----------------- //
  vector<complex<T>> sDelay(s.size());  // Our shifted spectrum.
  T phaseShift=(-2.0*M_PI*fShift/fs);   // Precompute phase shift.
  complex<T> delayMod(1.0,0.0);         // Initialize modulator to angle 0.
  for (size_t i=0; i<s.size(); ++i)     // For the existence of our signal..
  {
    sDelay[i]=s[i]*delayMod;            // Calculate this frequency bin.
    delayMod*=polar(1.0,phaseShift);    // Increment the phase.
  }                                     // Done delay shifting the signal.
  return sDelay;
}                                       // ---------- Shift ----------------- //

template<typename T>
// Method to perform a carrier sweep with a start and stop frequency about the 
// center frequency
vector<vector<complex<T>>> SpectralOps<T>::Sweep(
  const vector<complex<T>>&s,           // The input signal
  const double fStart,                  // The Start Frequency.
  const double fCenter,                 // The center frequency
  const double fStop,                   // The stop frequency
  const double step,                    // The number of bins to jump
  const double fs,                      // The sample rate of the signal
  const WindowType &w,                  // The window to apply to the signal
  const int wSiz,                       // The size of the window
  const float ovlap)                    // The overlap factor
{                                       // ---------- Sweep ----------------- //
    // -------------------------------- //
    // First we shift the input signal about the center frequency in 
    // our spectum down to 0 Hz. This normalizes our signal to simplify analysis.
    // -------------------------------- //
  vector<complex<T>> cSig=Shift(s,-fCenter, fs);
  vector<vector<complex<T>>> sMat; 
    // -------------------------------- //
    // Scan through different frequency bands relative to the center frequency.
    // -------------------------------- //
  for (double freq=fStart; freq<=fStop; freq+=step)
  {
    // -------------------------------- //
    // Having our signal at 0 Hz we observe the behaviour of our signal when 
    // shifted by various frequencies relative to the central moment. 
    // -------------------------------- //
    vector<complex<T>> sShift=Shift(cSig,freq,fs);
    // -------------------------------- //
    // Apply the window to the signal to reduce spectral leakage.
    // -------------------------------- //
    vector<complex<T>> S=OLAProcessor(sShift,w,wSiz,ovlap);
    // -------------------------------- //
    // Perform an FFT on the windowed signal to analyze the energy 
    // of the signal at this frequency offset.
    // -------------------------------- //
    vector<complex<T>> spectrum=FFT(sWindowed);
    // -------------------------------- //
    // Normalize the FFT result by a suitable factor (maybe wSiz?) This should
    // distribute the energy across the window.
    // -------------------------------- //
    for (size_T i=0;i<spectrum.size();++i)
      spectrum[i]/=static_cast<T>(wSiz);
    // -------------------------------- //
    // Next we recollect the resulting FFT spectrum for this frequency offset,
    // to obtain the full spectra that represents the signal's behaviour accross
    // the entire sweep range.
    // -------------------------------- //
    sMat.push_back(spectrum);
  }
    // -------------------------------- //
    // Having collected all of the signal's behaviour across the frequency 
    // sweep, return the full signal's spectra.
    // -------------------------------- //
  return sMat;
} 
#endif // SPECTRALOPS_H

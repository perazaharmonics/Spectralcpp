// ====================================================================================
// Description: A collection of common Spectral Transforms and manipulation methods.
//   used for various forms of Digital Signal Processing.
//
// Author: J. Enrique Peraza, JEP
// =======================================================================================

#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <complex>
#include <thread>
#include <future>
#include <chrono>
#include <stdexcept>
#include "DSPWindows.h"

namespace dsp::spectral
{
  using namespace std;
 
// Class template for spectral operations
template<typename T>
class WaveletOps
{
  // The wavelets we know.
  enum class WaveletType
  {
    Haar,                            // The Haar Wavelet type.
    Db1,                             // Daubechies wavelet type 1.
    Db6,                             // Daubechies wavelet type 6.
    Dym5,                            // Symlet type 5
    Sym8,                            // Symlet type 8.
    Coif5                            // Coiflet type 5 
  };
  enum class ThresholdType
  {
    Hard,                           // Hard thresholding.
    Soft                            // Soft thresholding.
  };
  public:
    // Constructors
    explicit WaveletOps(          // Constructor
      WaveletType wt=WaveletType::Haar, // Our default wavelet.
      size_t levels=1,            // Decomposition level.
      double threshold=0.0f,      // Denoising threshold.
      ThresholdType tt=ThresholdType::Hard)// The denoising type.
     : wType{wt},levels{levels},threshold{threshold},tType{tt} {}
    // ---------------------------- //
    // Denoise: Apply multi-level DWT denoising and reconstruct signal.
    // ---------------------------- //
    std::vector<double> Denoise(const std::vector<double>& signal)
    {                               // --------- Denoise --------------- //
      // -------------------------- //
      // 1. Zero-Pad the signal to make sure the operation is possible.
      // -------------------------- //   
      auto padded=pad_to_pow2(signal);// Zero-pad the signal.
      // -------------------------- //
      // 2. Now forward DWT + threshold.
      // -------------------------- //
      auto coeffs=dwt_multilevel(padded,selectForward(),
        this->levels,this->threshold,tTypeToString());
      // -------------------------- //
      // 3. Now perform the inverse DWT.
      // ------------------------- //
      auto recon=idwt_multilevel(coeffs,selectInverse());
      // ------------------------- //
      // 4. Remove the zero padding from the signal.
      // ------------------------- //
      return remove_padding(recon,signal.size());
    }                              // ---------- Denoise ----------- 
public:
  inline vector<double> pad_to_pow2(const vector<double>& signal) 
  {
    size_t original_length = signal.size();
    size_t padded_length = static_cast<size_t>(next_power_of_2(original_length));
    vector<double> padded_signal(padded_length);

    copy(signal.begin(), signal.end(), padded_signal.begin());
    fill(padded_signal.begin() + original_length, padded_signal.end(), 0);

    return padded_signal;
  }                                                     
/// Remove padding back to the original length
  inline vector<double> remove_padding(const vector<double>& signal, size_t original_length) 
  {
    return vector<double>(signal.begin(), signal.begin() + original_length);
  }
  
// Normalization.
inline vector<double> normalize_minmax(const vector<double>& data) 
{
    double min_val = *min_element(data.begin(), data.end());
    double max_val = *max_element(data.begin(), data.end());
    vector<double> normalized_data(data.size());

    transform(data.begin(), data.end(), normalized_data.begin(),
        [min_val, max_val](double x) { return (x - min_val) / (max_val - min_val); });

    return normalized_data;
}

inline vector<double> normalize_zscore(const vector<double>& data) 
{
    double mean_val = accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double std_val = sqrt(sq_sum / data.size() - mean_val * mean_val);
    vector<double> normalized_data(data.size());

    transform(data.begin(), data.end(), normalized_data.begin(),
        [mean_val, std_val](double x) { return (x - mean_val) / std_val; });

    return normalized_data;
}

inline vector<double> awgn(const vector<double>& signal, double desired_snr_db) 
{
    double signal_power = accumulate(signal.begin(), signal.end(), 0.0,
        [](double sum, double val) { return sum + val * val; }) / signal.size();
    double desired_noise_power = signal_power / pow(10, desired_snr_db / 10);
    vector<double> noise(signal.size());

    generate(noise.begin(), noise.end(), [desired_noise_power]() {
        return sqrt(desired_noise_power) * ((double)rand() / RAND_MAX * 2.0 - 1.0);
    });

    vector<double> noisy_signal(signal.size());
    transform(signal.begin(), signal.end(), noise.begin(), noisy_signal.begin(), plus<double>());

    return noisy_signal;
}
/// Multi‐level discrete wavelet transform with hard/soft thresholding
inline vector<pair<vector<double>, vector<double>>> 
  dwt_multilevel(vector<double>& signal, 
    function<pair<vector<double>, vector<double>>(const vector<double>&)> wavelet_func, 
    size_t levels, double threshold, const string& threshold_type = "hard") 
    {

      size_t n = signal.size();
      size_t n_pad = static_cast<size_t>(next_power_of_2(n));
      if (n_pad != n)
        signal.resize(n_pad, 0.0);
      vector<pair<vector<double>, vector<double>>> coeffs;
      vector<double> current_signal = signal;
      for (size_t i = 0; i < levels; ++i) 
      {
         auto [approx, detail] = wavelet_func(current_signal);
        // Apply thresholding to the detail coefficients
        if (threshold_type == "hard") 
          detail = hard_threshold(detail, threshold);
        else if (threshold_type == "soft")
          detail = soft_threshold(detail, threshold);
        coeffs.emplace_back(approx, detail);
        current_signal = approx;

        if (current_signal.size() < 2)
          break;
    }
    return coeffs;
}                                                                                                                         
/// Inverse multi‐level DWT
inline vector<double> 
idwt_multilevel(vector<pair<vector<double>, vector<double>>>& coeffs, function<vector<double>(const vector<double>&, const vector<double>&)> wavelet_func) 
{                                   // ------- idwt_multilevel ----
    vector<double> signal = coeffs[0].first;
    for (int i = coeffs.size() - 1; i >= 0; --i) 
    {
      auto [approx, detail] = coeffs[i];
      signal = wavelet_func(approx, detail);
    }

    return signal;
}

// Normalization algorithms.

/// Forward wavelet bases
inline pair<vector<double>, vector<double>> haar(const vector<double>& signal) 
{
    const vector<double> h = { 1 / sqrt(2), 1 / sqrt(2) };
    const vector<double> g = { 1 / sqrt(2), -1 / sqrt(2) };

    size_t n = signal.size() / 2;
    vector<double> approx(n), detail(n);

    for (size_t i = 0; i < n; ++i) {
        approx[i] = h[0] * signal[2 * i] + h[1] * signal[2 * i + 1];
        detail[i] = g[0] * signal[2 * i] + g[1] * signal[2 * i + 1];
    }

    return make_pair(approx, detail);
}                                     
inline pair<vector<double>, vector<double>> db1(const vector<double>& signal) 
{
    const vector<double> h = {
        (1 + sqrt(3)) / 4, (3 + sqrt(3)) / 4, (3 - sqrt(3)) / 4, (1 - sqrt(3)) / 4
    };
    const vector<double> g = { h[3], -h[2], h[1], -h[0] };

    size_t n = signal.size() / 2;
    vector<double> approx(n), detail(n);

    for (size_t i = 0; i < n; ++i) {
        approx[i] = 0;
        detail[i] = 0;
        for (size_t k = 0; k < h.size(); ++k) 
        {
            size_t index = 2 * i + k - h.size() / 2 + 1;
            if (index < signal.size()) {
                approx[i] += h[k] * signal[index];
                detail[i] += g[k] * signal[index];
            }
        }
    }

    return make_pair(approx, detail);
}                                     
inline pair<vector<double>, vector<double>> db6(const vector<double>& signal) 
{
    const vector<double> h = {
        -0.001077301085308,
        0.0047772575109455,
        0.0005538422011614,
        -0.031582039318486,
        0.027522865530305,
        0.097501605587322,
        -0.129766867567262,
        -0.226264693965440,
        0.315250351709198,
        0.751133908021095,
        0.494623890398453,
        0.111540743350109
    };
    const vector<double> g = {
        h[11], -h[10], h[9], -h[8], h[7], -h[6], h[5], -h[4], h[3], -h[2], h[1], -h[0]
    };

    size_t n = signal.size() / 2;
    vector<double> approx(n), detail(n);

    for (size_t i = 0; i < n; ++i) {
        approx[i] = 0;
        detail[i] = 0;
        for (size_t k = 0; k < h.size(); ++k) {
            size_t index = 2 * i + k - h.size() / 2 + 1;
            if (index < signal.size()) {
                approx[i] += h[k] * signal[index];
                detail[i] += g[k] * signal[index];
            }
        }
    }

    return make_pair(approx, detail);
}                                     
inline pair<vector<double>, vector<double>> sym5(const vector<double>& signal) 
{
    const vector<double> h = {
        0.027333068345078, 0.029519490925774, -0.039134249302383,
        0.199397533977394, 0.723407690402421, 0.633978963458212,
        0.016602105764522, -0.175328089908450, -0.021101834024759,
        0.019538882735287
    };
    const vector<double> g = {
        h[9], -h[8], h[7], -h[6], h[5], -h[4], h[3], -h[2], h[1], -h[0]
    };

    size_t n = signal.size() / 2;
    vector<double> approx(n), detail(n);

    for (size_t i = 0; i < n; ++i) {
        approx[i] = 0;
        detail[i] = 0;
        for (size_t k = 0; k < h.size(); ++k) {
            size_t index = 2 * i + k - h.size() / 2 + 1;
            if (index < signal.size()) {
                approx[i] += h[k] * signal[index];
                detail[i] += g[k] * signal[index];
            }
        }
    }

    return make_pair(approx, detail);
}
                                    
inline pair<vector<double>, vector<double>> sym8(const vector<double>& signal) 
{
    const vector<double> h = {
        -0.003382415951359, -0.000542132331635, 0.031695087811492,
        0.007607487325284, -0.143294238350809, -0.061273359067908,
        0.481359651258372, 0.777185751700523, 0.364441894835331,
        -0.051945838107709, -0.027219029917056, 0.049137179673476,
        0.003808752013890, -0.014952258336792, -0.000302920514551,
        0.001889950332900
    };
    const vector<double> g = {
        h[15], -h[14], h[13], -h[12], h[11], -h[10], h[9], -h[8],
        h[7], -h[6], h[5], -h[4], h[3], -h[2], h[1], -h[0]
    };

    size_t n = signal.size() / 2;
    vector<double> approx(n), detail(n);

    for (size_t i = 0; i < n; ++i) {
        approx[i] = 0;
        detail[i] = 0;
        for (size_t k = 0; k < h.size(); ++k) {
            size_t index = 2 * i + k - h.size() / 2 + 1;
            if (index < signal.size()) {
                approx[i] += h[k] * signal[index];
                detail[i] += g[k] * signal[index];
            }
        }
    }

    return make_pair(approx, detail);
}

inline pair<vector<double>, vector<double>> coif5(const vector<double>& signal) 
{
    const vector<double> h = {
        -0.000720549445364, -0.001823208870703, 0.005611434819394,
        0.023680171946334, -0.059434418646456, -0.076488599078311,
        0.417005184421393, 0.812723635445542, 0.386110066821162,
        -0.067372554721963, -0.041464936781959, 0.016387336463522
    };
    const vector<double> g = {
        h[11], -h[10], h[9], -h[8], h[7], -h[6], h[5], -h[4],
        h[3], -h[2], h[1], -h[0]
    };

    size_t n = signal.size() / 2;
    vector<double> approx(n), detail(n);

    for (size_t i = 0; i < n; ++i) {
        approx[i] = 0;
        detail[i] = 0;
        for (size_t k = 0; k < h.size(); ++k) {
            size_t index = 2 * i + k - h.size() / 2 + 1;
            if (index < signal.size()) {
                approx[i] += h[k] * signal[index];
                detail[i] += g[k] * signal[index];
            }
        }
    }

    return make_pair(approx, detail);
}
                                   
/// Inverse wavelet reconstruction
inline vector<double> inverse_haar(const vector<double>& approx, const vector<double>& detail) 
{
    const vector<double> h_inv = { 0.7071067811865476, 0.7071067811865476 };
    const vector<double> g_inv = { -0.7071067811865476, 0.7071067811865476 };

    vector<double> reconstructed_signal;
    for (size_t i = 0; i < approx.size(); ++i) {
        reconstructed_signal.push_back(approx[i] * h_inv[0] + detail[i] * g_inv[0]);
        reconstructed_signal.push_back(approx[i] * h_inv[1] + detail[i] * g_inv[1]);
    }

    return reconstructed_signal;
}
inline vector<double> inverse_db1(const vector<double>& approx, const vector<double>& detail)
{
    const vector<double> h_inv = {
        (1 + sqrt(3)) / 4, (3 + sqrt(3)) / 4, (3 - sqrt(3)) / 4, (1 - sqrt(3)) / 4
    };
    const vector<double> g_inv = { h_inv[3], -h_inv[2], h_inv[1], -h_inv[0] };

    vector<double> reconstructed_signal(2 * approx.size(), 0.0);
    for (size_t i = 0; i < approx.size(); ++i) {
        for (size_t k = 0; k < h_inv.size(); ++k) {
            size_t index = (2 * i + k - h_inv.size() / 2 + 1);
            if (index < reconstructed_signal.size()) {
                reconstructed_signal[index] += approx[i] * h_inv[k] + detail[i] * g_inv[k];
            }
        }
    }

    return reconstructed_signal;
}
inline vector<double> inverse_db6(const vector<double>& approx, const vector<double>& detail) 
{
    const vector<double> h_inv = {
        0.111540743350109, 0.494623890398453, 0.751133908021095,
        0.315250351709198, -0.226264693965440, -0.129766867567262,
        0.097501605587322, 0.027522865530305, -0.031582039318486,
        0.0005538422011614, 0.0047772575109455, -0.001077301085308
    };
    const vector<double> g_inv = {
        h_inv[11], -h_inv[10], h_inv[9], -h_inv[8], h_inv[7], -h_inv[6],
        h_inv[5], -h_inv[4], h_inv[3], -h_inv[2], h_inv[1], -h_inv[0]
    };

    vector<double> reconstructed_signal(2 * approx.size(), 0.0);
    for (size_t i = 0; i < approx.size(); ++i) {
        for (size_t k = 0; k < h_inv.size(); ++k) {
            size_t index = (2 * i + k - h_inv.size() / 2 + 1);
            if (index < reconstructed_signal.size()) {
                reconstructed_signal[index] += approx[i] * h_inv[k] + detail[i] * g_inv[k];
            }
        }
    }

    return reconstructed_signal;
}

inline vector<double> inverse_sym5(const vector<double>& approx, const vector<double>& detail) 
{
    const vector<double> h_inv = {
        0.019538882735287, -0.021101834024759, -0.175328089908450,
        0.016602105764522, 0.633978963458212, 0.723407690402421,
        0.199397533977394, -0.039134249302383, 0.029519490925774,
        0.027333068345078
    };
    const vector<double> g_inv = {
        h_inv[9], -h_inv[8], h_inv[7], -h_inv[6], h_inv[5], -h_inv[4],
        h_inv[3], -h_inv[2], h_inv[1], -h_inv[0]
    };

    vector<double> reconstructed_signal(2 * approx.size(), 0.0);
    for (size_t i = 0; i < approx.size(); ++i) {
        for (size_t k = 0; k < h_inv.size(); ++k) {
            size_t index = (2 * i + k - h_inv.size() / 2 + 1);
            if (index < reconstructed_signal.size()) {
                reconstructed_signal[index] += approx[i] * h_inv[k] + detail[i] * g_inv[k];
            }
        }
    }

    return reconstructed_signal;
}

inline vector<double> inverse_sym8(const vector<double>& approx, const vector<double>& detail) 
{
    const vector<double> h_inv = {
        0.001889950332900, -0.000302920514551, -0.014952258336792,
        0.003808752013890, 0.049137179673476, -0.027219029917056,
        -0.051945838107709, 0.364441894835331, 0.777185751700523,
        0.481359651258372, -0.061273359067908, -0.143294238350809,
        0.007607487325284, 0.031695087811492, -0.000542132331635,
        -0.003382415951359
    };
    const vector<double> g_inv = {
        h_inv[15], -h_inv[14], h_inv[13], -h_inv[12], h_inv[11], -h_inv[10],
        h_inv[9], -h_inv[8], h_inv[7], -h_inv[6], h_inv[5], -h_inv[4],
        h_inv[3], -h_inv[2], h_inv[1], -h_inv[0]
    };

    vector<double> reconstructed_signal(2 * approx.size(), 0.0);
    for (size_t i = 0; i < approx.size(); ++i) {
        for (size_t k = 0; k < h_inv.size(); ++k) {
            size_t index = (2 * i + k - h_inv.size() / 2 + 1);
            if (index < reconstructed_signal.size()) {
                reconstructed_signal[index] += approx[i] * h_inv[k] + detail[i] * g_inv[k];
            }
        }
    }

    return reconstructed_signal;
}

inline vector<double> inverse_coif5(const vector<double>& approx, const vector<double>& detail) 
{
    const vector<double> h_inv = {
        0.016387336463522, -0.041464936781959, -0.067372554721963,
        0.386110066821162, 0.812723635445542, 0.417005184421393,
        -0.076488599078311, -0.059434418646456, 0.023680171946334,
        0.005611434819394, -0.001823208870703, -0.000720549445364
    };
    const vector<double> g_inv = {
        h_inv[11], -h_inv[10], h_inv[9], -h_inv[8], h_inv[7], -h_inv[6],
        h_inv[5], -h_inv[4], h_inv[3], -h_inv[2], h_inv[1], -h_inv[0]
    };

    vector<double> reconstructed_signal(2 * approx.size(), 0.0);
    for (size_t i = 0; i < approx.size(); ++i) {
        for (size_t k = 0; k < h_inv.size(); ++k) {
            size_t index = (2 * i + k - h_inv.size() / 2 + 1);
            if (index < reconstructed_signal.size()) {
                reconstructed_signal[index] += approx[i] * h_inv[k] + detail[i] * g_inv[k];
            }
        }
    }

    return reconstructed_signal;
}

inline vector<double> hard_threshold(const vector<double>& detail, double threshold) 
{
    vector<double> result(detail.size());
    transform(detail.begin(), detail.end(), result.begin(), [threshold](double coeff) {
        return abs(coeff) < threshold ? 0.0 : coeff;
    });
    return result;
}

inline vector<double> soft_threshold(const vector<double>& detail, double threshold) 
{
    vector<double> result(detail.size());
    transform(detail.begin(), detail.end(), result.begin(), [threshold](double coeff) {
        return signbit(coeff) ? -max(0.0, abs(coeff) - threshold) : max(0.0, abs(coeff) - threshold);
    });
    return result;
}  
  
private:
  WaveletType wType{WaveletType::Haar};       // The wavelet of choice.
  ThresholdType tType{ThresholdType::Hard};     // The type of thresholding.
  size_t levels{1};              // The wavelet decomposition level.
  double threshold{0.0f};        // The threshold of when to cancel.
// Method to choose the correct forward wavelet function.
inline std::function<std::pair<std::vector<double>, std::vector<double>>(const std::vector<double>&)> selectForward(void) const
  {                              // -------- selectForward --------
    switch (wType)               // Act according to the type.
    {                            //
      case WaveletType::Haar: return haar;// We want Haar.
      case WaveletType::Db1:  return db1;// We want db1.
      case WaveletType::Db6: return db6;// We want db6.
      case WaveletType::Sym5: return sym5;// We want sym 5
      case WaveletType::Sym8: return sym8;// We want sym 8
      case WaveletType::Coif5: return coif5;// We want sym 5
    }                           // Done acting according to wlet typ.
    return haar;                // Return our default wavelet.
  }                             // -------- selectForward --------
// Chooses the correct inverse wavelet reconstruction
  inline std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)>
  selectInverse(void) const   // Select the correct reconstruct wave.
  {                           // -------- selectInverse --------
    switch(wType)             // Act according to the wave type
    {                         // Select the recon wavelet.
      case WaveletType::Haar: return inverse_haar;// We want Haar.
      case WaveletType::Db1:  return inverse_db1;// We want db1.
      case WaveletType::Db6: return inverse_db6;// We want db6.
      case WaveletType::Sym5: return inverse_sym5;// We want sym 5
      case WaveletType::Sym8: return inverse_sym8;// We want sym 8
      case WaveletType::Coif5: return inverse_coif5;// We want sym 5 
    }                          // Done acting according to wtype.
    return inverse_haar;       // Return inverse of default fwd wave.
  }                            // -------- selectInverse --------
  // Convert enum to the string the denoiser expects.
  inline std::string tTypeToString(void) const
  {
    return (this->tType==ThresholdType::Hard?"hard":"soft");
  }
protected:
/// Pad a signal up to the next power of two
  inline double next_power_of_2(double x) 
  {
    return x == 0 ? 1 : pow(2, ceil(log2(x)));
  }
  
};
// ----------------------------------------------------------------------------
// The discrete Cosine Transform (DCT) part of FCWTransforms.h
// Discrete Cosine / Modified Discrete Cosine Transform. 
// Supports DCT-I, DCT-II, DCT-III, and DCT-IV (1-D, real-valued).
// Supports MDCT/IMDCT (Princen-Bradley sine window pair).
// Uses the existing SpectralOps<T>::FFTStride/IFFTStride for O(N*log(N))
// No heap allocations inside the hot path (caller passes working buffers).
// --------------------------------------------------------------------------- 
template<typename T=float>
class DCT
{
  public:
    enum class Type { I, II, III, IV,MDCT,IMDCT };

    // ----------- Forward 1-D transforms ----------------------------------- 
    static vector<T> Transform(
      const vector<T>& x,               // The input signal to transform.
      Type t,                           // The transform type.
      SpectralOps<T>& engine)           // Our spectral engine.
    {                                   // ----------- Transform -------------
      switch (t)                        // Act according to transform type...
      {                                 // 
        case Type::I:     return DCTI(x,engine); // DCT-I
        case Type::II:    return DCTII(x,engine); // DCT-II
        case Type::III:   return DCTIII(x,engine); // DCT-III
        case Type::IV:    return DCTIV(x,engine); // DCT-IV
        case Type::MDCT:  return MDCT(x,engine); // MDCT
        case Type::IMDCT: return IMDCT(x,engine); // IMDCT
        default:          return DCTIV(x,engine); // Default to DCT-II
      }                                 // Done dispatching transform          
    }                                   // ----------- Transform -------------
    static vector<T> Forward(
      const vector<T>& timeBlock,       // The time-domain block to transform.
      const Window<T>& analysisW,       // The analysis window to apply.
       SpectralOps<T>& engine)          // Our spectral engine.
    {                                   // ----------- Forward --------------- //
      return MDCT(timeBlock,analysisW,engine); // Use MDCT for forward transform.
    }                                   // ----------- Forward --------------- //
    // ----------- Inverse 1-D transforms -----------------------------------
    static vector<T> Inverse(
      const vector<T>& X,               // The frequency-domain block to transform.
      Type t,                           // The transform type.
    SpectralOps<T>& engine)             // Our spectral engine.
    {                                   // ----------- Inverse ----------------
      switch(t)                         // Dispatch according to type.
      {                                 //
        case Type::I:     return DCTI(X,engine,/*inverse=*/true); // DCT-I
        case Type::II:    return DCTIII(X,engine); // -1 dual
        case Type::III:   return DCTII(X,engine); // -1 dual
        case Type::IV:    return DCTIV(X,engine); // Self-inverse
        default:          return DCTIV(X,engine); // Default to DCT-II
      }                                 // Done dispatching inverse transform
    }                                   // ----------- Inverse ----------------
    static vector<T> Inverse(
      const vector<T>* coeffs,          // The frequency-domain coefficients to transform.
      const Window<T>& synthesisW,      // The synthesis window to apply.
    SpectralOps<T>& engine)             // Our spectral engine.      
    {                                   // ----------- Inverse ----------------
      return IMDCT(coeffs,synthesisW,engine); // Use IMDCT for inverse transform.
    }                                   // ----------- Inverse ----------------
    // ---------------- Discrete Cosine Transforms ----------------
    // All four classical DCTs are computed through a single length-2N complex FFT.
    // following the well-known even/odd embeddings (see Britanak & Rao, ch. 6).
    // We use the existing SpectralOps<T>::FFTStride/IFFTStride for O(N*log(N))
    // So no extra radix-2 code is duplicated here.
    // 
    static vector<T> DCTII(
      const vector<T>& x,               // The input signal to transform.
      SpectralOps<T>& engine,           // Our spectral engine.
      const bool inverse=false)         // If it is the inverse.           
    {                                   // ----------- DCT-II ----------------
      const size_t N=x.size();          // Get the size of the input signal.
      vector<std::complex<T>> w(2*N);     // Create a complex vector of size 2*N.
      // Even-symmetrical extensions (x, reversed(x))
      for (size_t n=0;n<N;++n)          // For all samples.
      {                                 // Set the even-symmetrical extensions.
        w[n]={x[n],0};                  // Fill the first N points with x[n].
        w[2*N-1-n]={x[n],0};            // Fill the last N points with x[n].
      }                                 // Done symmetrical extensions.
      // Compute the FFT of the complex vector w.
      auto W=engine.FFTStride(w);       // Compute the FFT of the complex vector w.
      // Post-twiddle & packing.
      vector<T> X(N);                   // Create a vector of size N for the output.
      const T scale=inverse?2./N:1.;    // Orthonormal pair (II <--> III).
      const T factor=M_PI/(2.*N);       // The factor to apply to the output.
      for (size_t k=0;k<N;++k)          // For all frequency bins...
      {                                 // Apply the post-twiddle and scale/pack.
        std::complex<T> c=std::exp(std::complex<T>(0,-factor*k));
        X[k]=2*(W[k]*c).real()*scale;   // Apply the post-twiddle and scale.
      }                                 // Done post-twiddle and packing.
      return X;                         // Return the DCT-II coefficients.
    }                                   // ----------- DCT-II ----------------
    static std::vector<T> DCTIII(
      const std::vector<T>& x,          // The input signal to transform.
      SpectralOps<T>& engine)           // Our spectral engine.
    {                                   // ----------- DCT-III ---------------
      // DCT-III is the inverse of DCT-II. Call DCTII with the inverse flag
      return DCTII(x,engine,true);      // Call DCT-II with the inverse flag.
    }                                   // ----------- DCT-III ---------------
    static vector<T> DCTI(
      const std::vector<T>& x,          // The input signal to transform.
      SpectralOps<T>& engine,           // Our spectral engine.
      const bool inverse=false)         // If it is the inverse.
    {                                   // ----------- DCT-I -----------------
      // Even-even extensions --> length-2(N-1) FFT.
      const size_t N=x.size();          // Get the size of the input signal.
      if (N<2) return x;                // If the input signal is too short, return it.
      vector<std::complex<T>> w(2*(N-1)); // Create a complex vector of size 2*(N-1).
      // Fill the complex vector with the even-even extensions.
      for (size_t n=0;n<N-1;++n)        // For all samples.
      {                                 // Set the even-even extensions.
        w[n]={x[n],0};                  // Fill the first N-1 points with x[n].
        w[2*(N-1)-1-n]={x[N-2-n],0};    // Fill the last N-1 points with x[N-2-n].
      }                                 // Done even-even extensions.
      w[N-1]={x.back(),0};              // Fill the middle point with x[N-1].
      auto W=engine.FFTStride(w);       // Compute the FFT of the complex vector w.
      vector<T> X(N);                   // Create a vector of size N for the output.
      const T scale=inverse?1./(N-1):1.; // Orthonormal pair (I <--> IV).
      for (size_t k=0;k<N;++k)          // For all frequency bins...
        X[k]=W[k].real()*scale;         // Apply the post-twiddle and scale.
      return X;                         // Return the DCT-I coefficients.
    }                                   // ----------- DCT-I -----------------
    static vector<T> DCTIV(
      const vector<T>& x,               // The input signal to transform.
      SpectralOps<T>& engine)           // Our spectral engine.
    {                                   // ----------- DCT-IV ----------------
      const size_t N=x.size();          // Get the size of the input signal.
      vector<std::complex<T>> w(2*N); // Create a complex vector of size 2*N.
      // Fill the complex vector with the even-even extensions.
      for (size_t n=0;n<N;++n)          // For all samples...
      {
        w[n]={x[n],0};                  // Fill the first N points with x[n].
        w[2*N-1-n]={-x[n],0};           // Fill the last N points with -x[n].
      }                                 // Done symmetrical extensions.
      auto W=engine.FFTStride(w);       // Compute the FFT of the complex vector w.
      vector<T> X(N);                   // Create a vector of size N for the output.
      const T factor=M_PI/(4.*N);       // The factor to apply to the output.
      for (size_t k=0;k<N;++k)          // For all frequency bins....
      {                                 
        std::complex<T> c=std::exp(std::complex<T>(0,-factor*(2*k+1)));
        X[k]=2*(W[k].real()*c.real()-W[k].imag()*c.imag());// Apply the post-twiddle and scale.
      }                                 // Done post-twiddle and packing.
      return X;                         // Self-inverse (orthogonal, no extra scaling).
    }                                   // ----------- DCT-IV ----------------
  // --------------------- MDCT/IMDCT ---------------------
  // MDCT length is N (produces N/2 coefficients from a 2N-sample block).
  // IMDCT returns a 2N block that the caller overlap-adds by N samples.
  // -----------------------------------------------------------
  static vector<T> MDCT(
    const vector<T>& timeBlock,         // The time-domain block to transform.
    const Window<T>& win,               // The window to apply before MDCT.
    SpectralOps<T>& engine)             // Our spectral engine.
  {                                     // ------------ MDCT ---------------
    // Preconditions:
    const size_t n2=timeBlock.size();   // Get the size of the input block.
    if (n2%2)                           // If a remainder appears after division by 2....
      return vector<T>();               // Return an empty vector.
    const size_t N=n2/2;                // N is the number of MDCT coefficients.
    // -------------------------------- //
    // Windowing + pre-twiddle (+ phase shift n+0.5).
    // -------------------------------- //
    vector<std::complex<T>> xwin(n2);   // Create a complex vector of size n2.
    for (size_t n=0;n<n2;++n)           // Fpr all samples...
      xwin[n]={timeBlock[n]*win[n],0}; // Apply the window to the time block.
    // -------------------------------- //
    // Rearrange into even-odd groups for FFT of length 2N.
    // -------------------------------- //
    vector<std::complex<T>> v(n2);      //  Create a complex vector of size n2.
    for (size_t n=0;n<N;++n)            // For all samples....
    {                                   // Fill the complex vector with the even-odd groups.
      v[n]=xwin[n]+xwin[n2-1-n];        // Even part: x[n]+w[n2-1-n].
      v[n+N]=(xwin[n]-xwin[n2-1-n])*std::complex<T>{0,-1};     // Odd part: x[n]-w[n2-1-n].
    }                                   //
    auto V=engine.FFTStride(v);         // Compute the FFT of the complex vector v.
    // Post-twiddle: take real part, multiply by exp(-j(pi*k+0.5)/N)                                                
    vector<T> X(N);                    // Create a vector of size N for the output.
    const T factor=M_PI/(2.*N);         // The factor to apply to the output.
    for (size_t k=0;k<N;++k)            // For all frequency bins...
    {                                   // Apply the post-twiddle and scale.
      std::complex<T> c=std::exp(std::complex<T>(0,-factor*(k+0.5)));
      X[k]=(V[k]*c).real();             // Take the real part and apply the post-twiddle.
    }                                   // Done post-twiddle and packing.
    return X;                           // Return the MDCT coefficients.
  }                                     // ------------ MDCT ---------------
  static vector<T> IMDCT(
    const vector<T>& X,                 // The frequency-domain block to transform.
    const Window<T>& win,               // The window to apply after IMDCT.
    SpectralOps<T>& engine)             // Our spectral engine.
  {                                     // ------------ IMDCT ---------------
    const size_t N=X.size();            // Get the size of the input block.
    const size_t n2=N*2;                //
    // Pre-twiddle
    vector<std::complex<T>> V(n2);      // Create a complex vector of size n2.
    const T factor=M_PI/(2.*N);         // The factor to apply to the output.
    for (size_t k=0;k<N;++k)            // For all frequency bins....
    {
      std::complex<T> c=std::exp(std::complex<T>(0,factor*(k+0.5)));
      V[k]=X[k]*c;                      // Apply the pre-twiddle.
      V[k+N]=-X[k]*c;                   // odd symmetry: V[k+N]=-X[k]*c.
    }                                   // 
    // IFFT but actua;;y FFT of size 2N, then scale/flip.
    auto v=engine.IFFTStride(V);        // Compute the IFFT of the complex vector V.
    // Post-twiddle + windowing.
    vector<T> y(n2);                     // Output signal.
    for (size_t n=0;n<n2;++n)           // For all samples...
      y[n]=2*(v[(n+N/2)%n2].real())*win[n]; // Apply the post-twiddle and windowing.
    return y;           // Caller must perform 50% OLA with previous frame (block).
  }                                     // ------------ IMDCT ---------------
};  



// The Fourier part of FCWTransforms.h
template<typename T>
class SpectralOps
{
    using WindowType = typename Window<T>::WindowType; // Alias for WindowType
public:
    // Constructors


    // Accessors
    const inline vector<T> GetSignal (void) const { return signal; }
    inline void SetSignal (const vector<complex<T>>& s) {signal=s;}
    const inline int GetSamples (void) const { return length; }
    inline void SetSamples (const int N) {length=N;}
    const inline double GetSampleRate (void) const { return sRate; }
    inline void SetSampleRate (const double fs) {sRate=fs;}
    const inline vector<complex<T>> GetTwiddles (void) const { return twiddles; }
    inline void SetSubCarrier(const vector<complex<T>> &s) { subCarrier = s; }
    const inline vector<complex<T>> GetSubCarrier (void) { return subCarrier; } 
    
// ------------------------------------ //
// Constructors and Destructors
// ------------------------------------ //

SpectralOps(void)
{
  this->length = 0; // Set the length to 0.
  this->sRate = 0.0; // Set the sample rate to 0.
  this->window = WindowType::Hanning; // Default window type is Hamming.
  this->windowSize = 1024; // Default window size is 1024.
  this->overlap = 0.5; // Default overlap is 50%.
  this->signal.clear(); // Clear the signal vector.
  this->twiddles.clear(); // Clear the twiddles vector.
  this->subCarrier.clear(); // Clear the subcarrier vector.

}

SpectralOps(const vector<T> &s, const WindowType &w, const int windowSize) : signal{s}, length{0}, sRate{0.0}, window{w}, windowSize{windowSize}, overlap{0.5}
{
    // Generate the appropriate window for the given window size
    this->window = w.SetWindowType(window, windowSize);
}


SpectralOps(const vector<T> &s, const WindowType &w, const int windowSize, const float overlap) : signal{s}, length{0}, sRate{0.0}, window{w}, windowSize{windowSize}, overlap{overlap}
{
    // Generate the appropriate window for the given window size
    this->window=w.SetWindowType(window, windowSize);
}


~SpectralOps(void)
{
    signal.clear();
    twiddles.clear();
}


// ======================== Utility Methods ================================= //
// Utility Methods to precompute operations needed for Spectral Manipulations.
// ========================================================================== //

inline vector<complex<T>> TwiddleFactor(int N)
{
    if (twiddles.size() != N / 2)                        // Did we precompute N/2 twiddles before?
    {                                                    // No, so we..
        twiddles.resize(N / 2);                          // Resize the twiddles factor vector.
        for (int i = 0; i < N / 2; ++i)                  //  loop for the N/2 points and
            twiddles[i] = polar(1.0, -2 * M_PI * i / N); //  compute the twiddles factors.
    }
    return twiddles;
}
// Get the smallest power of 2 that is greater than or equal to N
// that can hold the input sequence for the Cooley-Tukey FFT,
// which splits the input sequence into even and odd halves.

inline int UpperLog2(const int N)
{
    for (int i = 0; i < 30; ++i) // For the first 30 powers of 2
    {                            // Compute the power of 2 as 2^i
      const int mask = 1 << i;   // Compute the value of 2^i
      if (mask >= N)             // If the power of 2 is >= N
        return i;                // Return the smallest power of 2 (i).
    }                            //
    return 30;                   // Else return 30 as the upper bound.
}


inline vector<int>ToInt(const vector<complex<T>> &s)
{
    vector<int> sInt(s.size());
    for (size_t i = 0; i < s.size(); ++i)
        sInt[i] = static_cast<int>(s[i].real());
    return sInt;
}


inline vector<double> ToReal(const vector<complex<T>> &s)
{
    vector<double> sReal(s.size());
    for (size_t i = 0; i < s.size(); ++i)
        sReal[i] = s[i].real();
    return sReal;
}
// Determine the amount of frequency bins to analyze per second of data.

inline vector<T> SetRBW(double rbw, double fs)
{
  const int wSiz=static_cast<int>(fs/rbw);
  // Window is assumed to have been defined by the caller before calling this
  // method.
  return GetWindow(window,wSiz);
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
/* @param twiddles: Vector of precomputed twiddles factors.
/* @param rot: The current stage of the FFT, iteration indicator. Starts at 0 for the first stage.
/* @param nBits: log2(N) where N is the length of the signal, total number of FFT stages.
/* Reference: https://github.com/AndaOuyang/FFT/blob/main/fft.cpp
*/

inline void ForwardButterfly(vector<T> &last, vector<T> &curr, const vector<T> &twiddles, const int rot, const int nBits)
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
    // Multiply the odd element by the twiddles factor for this butterfly group.  
    // The twiddles factor is a complex number that rotates the odd index.
    // and introduces the phase shift needed for the FFT. 
    // ------------------------------------- //   
        last[oddNdx] *= twiddles[j * phases];// Multiply the odd index by the twiddles factor.
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
template<typename U>
inline void ForwardButterfly(
  vector<std::complex<U>>& last,
  vector<std::complex<U>>& curr,
  const vector<std::complex<U>>& twiddles,
  int rot,int nBits)
{
  if (rot==nBits) return;
  const int sect=1<<(rot+1);            // Size of the butterfly section.
  const int phases=last.size()/sect;    // Number of sections the signal is divided into.
  for (int i=0;i<phases;++i)            // For every phase in the FFT
  {                                      // Perform the butterfly operation.
    const int base=i*sect;               // Base index for the current phase.
    for (int j=0;j<sect/2;++j)           // For every butterfly group in the structure.
    {
      const int evenNdx=base+j;          // Even index in the butterfly group.
      const int oddNdx=base+sect/2+j;    // Odd index in the butterfly group.
      last[oddNdx]*=twiddles[j*phases];  // Multiply the odd index by the twiddles factor.
      curr[evenNdx]=last[evenNdx]+last[oddNdx]; // Compute the even index.
      curr[oddNdx]=last[evenNdx]-last[oddNdx];  // Compute the odd index.
    }                                    // Done with all butterfly groups.
  }                                      // Done with all phases.
  ForwardButterfly(curr, last, twiddles, rot + 1, nBits); // Recurse to the next stage.
}
// Bit reversal permutation for the Cooley-Tukey FFT algorithm.

inline void  BitReversal(vector<T> &s, const int nBits)
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
// Overloaded for complex vectors.
inline void  BitReversal(vector<std::complex<T>> &s, const int nBits)
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
    
    inline std::pair<std::vector<T>, T>
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
    //    at fftsize uniformly spaced frequencies [0, 2π).  We solve AR(p) via
    //    Levinson‐Durbin, then evaluate
    //      H(ω) = σ² / |1 + a[1] e^{-jω} + … + a[p] e^{-j p ω} |²
    //    at Nfft points, returning a vector<complex<T>> of length Nfft
    //    (you can take real(H) or abs(H)² as your PSD). 
    // ------------------------------------------------------------------------
    
    inline std::vector<std::complex<T>>
    AR_PSD(const std::vector<T>& r, int order, int fftsize) const
    {
        if (order < 1 || (int)r.size() < order+1) {
            throw std::invalid_argument{"AR_PSD: order must be ≥1 and r.size() ≥ order+1"};
        }
        // 1) run Levinson‐Durbin on r[0..order]
        auto [a, sigma2] = LevinsonDurbin(r, order);
        // a = vector length p, contains a[1],…a[p], and sigma2 = error at final stage

        // 2) build PSD at fftsize freq bins
        std::vector<std::complex<T>> psd(fftsize);
        const T normFactor = T{2} * M_PI / static_cast<T>(fftsize);
        for (int k = 0; k < fftsize; ++k) {
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

inline std::pair<vector<complex<T>>,vector<vector<complex<T>>>> FFTStrideEig(const vector<complex<T>> &s)
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
  BitReversal(last, NBits);   // Perform bit reversal permutation.
  // ---------------------------------- //
  // Perform the FFT butterfly operation for the Cooley-Tukey FFT.
  // This computes the FFT in-place using the last and current buffers.
  // This is where the Cooley-Tukey FFT algorithm takes place.
  // ---------------------------------- //
  ForwardButterfly(last, curr, twiddles, 0, NBits); // Perform the FFT butterfly.
  // ---------------------------------- //
  // Here we compute the Fourier matrix and index the eigenvectors.
  // Return the computed FFT spectrum and the eigenvectors of the FFT matrix.
  // The eigenvectors are the Fourier basis vectors, which are the twiddles factors.
  // ---------------------------------- //
  vector<vector<complex<T>>> eigvecs(N, vector<complex<T>>(N)); // Create a matrix for the eigenvectors.
  const T invsqrt=static_cast<T>(1)/sqrt(static_cast<T>(N)); // Inverse square root of N for normalization.
  for (int ell=0;ell<N;ell++)           // For each row...
  {                                     // Compute the eigenvector for the row.
  // The row index ell corresponds to the frequency bin in the FFT.
      for (int k=0;k<N;k++)             // For each col in the eigenvector matrix...
    {                                   // Compute the Fourier matrix.
      long double angle=-2.0L*M_PI*(static_cast<long double>(ell))*(static_cast<long double>(k))/(static_cast<long double>(N));
      eigvecs[ell][k]=complex<T>(std::cos(angle),std::sin(angle))*invsqrt; // Compute the k-th eigenvector.
    }                                   // End of the loop.
  }                                     // End of the loop.
  return {curr, eigvecs};               // Return the computed FFT spectrum and the eigenvectors.
}


inline vector<complex<T>> FFTStride (const vector<complex<T>> &s)
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
    // Precompute the twiddles factors for the FFT.
    // The twiddles factors are used to rotate the signal in the FFT.
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
    BitReversal(last, nBits);   // Perform bit reversal permutation.
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
// The IFFT can be computed using the FFT with flipped order of the 
// frequency bins. That is, the complex conjugate of the input signal.
//   and thus the twiddles factors.
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

inline vector<complex<T>>  IFFTStride (const vector<complex<T>>& s)
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
// 3. Now we precompute the twiddles factors, this saves A LOT of time:
//      TwiddleFactor(k) = exp( -j * (2*pi/n)*k)
// 4. Having the twiddles Factors we compute the FFT butterfly to obtain the full
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

inline vector<complex<T>>  FFT(const vector<complex<T>>& s)
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
    // Precompute the twiddles factors.
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
    t=twiddles[k]*oddF[k];               // Scale and get this freq bin.
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
// The IFFT is nothing more than the FFT multiplied by 1/N, and with a twiddles
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

inline vector<complex<T>>  IFFT (const vector<complex<T>> &s)
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
  vector<complex<T>> sig(N);            // Signal buffer.
  for (int i=0;i<N;++i)                 // For all samples in reversed spectrum
    sig[i]=conj(S[i])/static_cast<T>(N);// Reverse and normalize.
  return sig;                          // Return time-domain signal.  
} 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Convolution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Method to perform the convolution of two signals. Typically in our context
// the convolution is done between an input signal s(n) of length N  and
// filter of length N h(n).
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

inline vector<complex<T>> Convolution( // Obtain the product of two signals.
    const vector<complex<T>> &s,        // Our input signal.
    const vector<complex<T>> &h)        // Our filter.
{                                       // ---------- Convolution ----------- //
    const int n = s.size();             // Length of the input signal.
    const int m = h.size();              // Length of the filter.
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

// Multiply two-length-N spectra element by element.

inline vector<complex<T>> PointwiseMul(const std::vector<std::complex<T>>& A,const std::vector<std::complex<T>>& B) const
{                                       // ---------- PointwiseMul ----------- //
  assert(A.size()==B.size());           // Ensure both spectra are of the same length.
  std::vector<std::complex<T>> C(A.size()); // Create a vector to hold the result.
  for (size_t i=0;i<A.size();++i)       // For each spectral element...
    C[i]=A[i]*B[i];                     // Multiply the two spectra element by element.
  return C;                             // Return the resulting spectrum.
}                                       // ---------- PointwiseMul ----------- //

// Short-Time Fourier Transform

inline vector<vector<complex<T>>> STFT(const vector<complex<T>> &s, 
  const WindowType &w, 
  int wSiz, 
  const float overlap)
{
  int step = wSiz * (1 - overlap / 100.0);// step size based on the ovlap %
  int nSegs = (s.size() - wSiz + step) / step;// # segs for the signal.
  vector<vector<complex<T>>> sMat(nSegs, vector<complex<T>>(wSiz));
  vector<T> window = this->window->GetWindow();// Window to be applied.
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
      sMat[i] = FFT(seg);               // Compute the FFT of the windowed segment.
  }                                     // Done with all segments.            
  return sMat;                          // The windowed spectrum.
}

// Inverse Short-Time-Fourier-Transform Method.
inline vector<complex<T>> ISTFT(
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
    vector<T> window = this->window.GenerateWindow(w, wSiz); // Get the window to be applied.
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

inline vector<complex<T>> OLAProcessor(
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

inline vector<complex<T>> OLAProcessor(
    const vector<complex<T>> &s, // The input signal.
    const WindowType &w,         // The window used.
    const int wSiz,              // The size of the window.
    const float ovlap)           // The percentage of overlap
{
    vector<vector<complex<T>>> sMat = STFT(s, w, wSiz, ovlap);// STFT of input signal.
    vector<complex<T>> H(wSiz,complex<T>(0.0,0.0));// Dummy Impulse filter.
    H=FFT(H);                           // FFT of the FIR filter.
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

inline vector<T> WelchPSD(
 const vector<T> &s,                    // The signal to process.
 const WindowType& w,                   // The window to apply to the signal.
 const int wSiz,                        // The size of the window.
 const float ovlap,                     // Overlap percentage (50% typical)
 const int fftsiz)                      // The size of the FFT.
{
    // -------------------------------- //
    // Compute the STFT of the signal.
    // -------------------------------- //
  vector<vector<complex<T>>> stftMat=STFT(s,w,wSiz,ovlap);
    // -------------------------------- //
    // Determine the scale factor of the window.
    // -------------------------------- //
  vector<T> pxxAvg(length,T(0));             // Initialize PSD Buffer
  const double winScl=pow(norm(this->window.GenerateWindow(w,wSiz),2),2); // Get Scale Factor.
    // -------------------------------- //
    // Now we accumulate the PSD for each segment in the STFT matrix.
    // -------------------------------- //
  for (int i=0; i<stftMat.size();++i)   // For each fft segment.
  {                                     // Where each row in the STFT matrix is
    const vector<complex<T>>& fftSegment=stftMat[i];//  an FFT segment
    // -------------------------------- //
    // Compute the Power Spectal Density
    // -------------------------------- //
    for (int j=0;j<fftsiz;++j)          // For every sample to process by the FFT
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
  for (int i=0; i<fftsiz;++i)           // For every sample...
    pxxAvg[i]/=(2*M_PI);                //
    // -------------------------------- //
    // Ensure the total energy of the sigal is conserved in the sprectrm.
    // Parseval's Theorem: SUM[n to N-1] x[n]^2 = (1/N) SUM[n to N-1] X[k]^2
    // -------------------------------- //
  pxxAvg[0]/=2;                         // Avergave freq bin 1 (DC component).
  for (int i=0;i<fftsiz;++i)            // For ever sample count the energy in
    pxxAvg[i]*=2;                       //  positive and negative halves
    // -------------------------------- //
    // Return the first half of the PSD (our FFT is symmetric) and we already
    // have recollected all the power.
    // -------------------------------- //
  return vector<T>(pxxAvg.begin(),pxxAvg.being()+fftsiz/2+1);

}
// Method to perform a frequency shift of the center frequency by a const amount

inline vector<complex<T>> Shift(  // Shift the signal in frequency domain.
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
// Method to perform a carrier sweep with a start and stop frequency about the 
// center frequency

inline vector<vector<complex<T>>> Sweep(
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
    vector<complex<T>> spectrum=FFT(S);
    // -------------------------------- //
    // Normalize the FFT result by a suitable factor (maybe wSiz?) This should
    // distribute the energy across the window.
    // -------------------------------- //
    for (size_t i=0;i<spectrum.size();++i)
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

private:
    vector<T> signal;                        // The signal to process.
    vector<T> subCarrier;                    // 
    WindowType window;                       // The window to apply to the signal.
    int windowSize=0;                    // The size of the window.
    float overlap;                     // The overlap factor.
    vector<complex<T>> twiddles;             // Precomputed twiddles factors.
    double sRate;                      // The rate at which we sampled the RF.
    int length;                        // The length of the signal.
};


}

#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <valarray>
#include <cstdint>

namespace dsp::spectral
{
// Missing Theory
using namespace std;
template<typename T>
class SpectralOps;
template<typename T>
class Window
{
public:
    enum class WindowType
    {
        Hanning,
        Hamming,
        BlackmanHarris,
        ExactBlackman,
        Blackman,
        FlatTop,
        FourTermBHarris,
        SevenTermBHarris,
        LowSideLobe,
        Rectangular,
        Tukey,
        Bartlett,
        Gaussian,
        Kaiser,
    };

    Window(void)
    {
      this->windowsize = 1024; // Default window size.
      this->window = WindowType::Rectangular; // Default window type.
      this->data=Rectangular(windowsize); // Initialize the window data with a rectangular window.
      this->data.reserve(windowsize); // Reserve space for the window data.
      this->data.resize(windowsize); // Resize the data vector to the window size.
      
    }
    ~Window(void)=default;
    // Access the elements of the window using the [] operator.
    T operator[](const size_t idx) const
    {
        return data[idx]; // Return the element at index idx.
    }
    T& operator[](size_t idx) {return data[idx];}
    auto begin(void) noexcept { return data.begin(); }
    auto end(void) noexcept { return data.end(); }
    auto begin(void) const noexcept { return data.cbegin(); }
    auto end(void) const noexcept { return data.cend(); }
    // Accessors
   inline  const Window<T> GetWindow (void) const { return window; }
   const inline size_t GetWindowsize (void) const { return windowsize; }
   inline vector<T> GetDefaultWindow (void) { return Rectangular(windowsize);}  
   inline void SetWindowsize (const size_t wsiz) {windowsize=wsiz;}
   inline WindowType GetWindowType (void) const { return window; }   
   inline void SetAlpha (const T a) { alpha=a;GenerateWindow(window, windowsize, a); } // Set the alpha value for Tukey window.
   inline T GetAlpha (void) const { return alpha; } // Get the alpha value for Tukey window.
   inline void SetSigma(const T s) { sigma=s; GenerateWindow(window, windowsize, alpha, s); } // Set the sigma value for Gaussian window.
   inline vector<T> GetData (void) const { return data; } // Get the window data.
   inline size_t Size(void) const { return data.size(); } // Get the size of the window data.
   inline void Clear(void) { data.clear(); } // Clear the window data.
   inline void Resize(const size_t N) { data.resize(N); } // Resize the window data to N elements.
   inline void Reserve(const size_t N) { data.reserve(N); } // Reserve space for N elements in the window data.
    // ------------------------------------ //
    // Constructors and Destructors
    // ------------------------------------ //
    Window(const size_t N)
        : windowsize(N), window(WindowType::Rectangular) {}
    Window(const WindowType &w, const size_t N) 
    {
    SetWindowType(w,N);
    }

    inline void SetWindowType(const WindowType &w, const size_t N)
    {
    window=w;                            // Store the type
    windowsize=N;                        // This long.
    GenerateWindow(w, N);               // Generate the window.
    }
    inline void SetWindowType(const WindowType &w, const size_t N, const T alpha)
    {
        window=w;                        // Store the type
        windowsize=N;                    // This long.
        this->alpha=alpha;              // Set the alpha value for Tukey window.
        GenerateWindow(w, N);           // Generate the window.
    }
    inline vector<T> GenerateWindow(const WindowType& w, const size_t N,const T alpha=0.5, const T sigma=0.4)
    {
        switch (w)                          // Set windows according to window type.
        {
            case WindowType::Hanning:         data=Hanning(N);break;
            case WindowType::Hamming:         data=Hamming(N);break;
            case WindowType::BlackmanHarris:  data=BlackmanHarris(N);break;
            case WindowType::ExactBlackman:   data=ExactBlackman(N);break;
            case WindowType::Blackman:        data=Blackman(N);break;
            case WindowType::FlatTop:         data=FlatTop(N);break;
            case WindowType::FourTermBHarris: data=FourTermBHarris(N);break;
            case WindowType::SevenTermBHarris:data=SevenTermBHarris(N);break;
            case WindowType::LowSideLobe:     data=LowSideLobe(N);break;
            case WindowType::Rectangular:     data=Rectangular(N);break;
            case WindowType::Tukey:           data=Tukey(N, alpha);break;
            case WindowType::Bartlett:        data=Bartlett(N);break;            
            case WindowType::Gaussian:        data=Gaussian(N,sigma);break;
            case WindowType::Kaiser:          data=Kaiser(N, alpha);break;
            default:                          data=Rectangular(N);
        }
        return data;                          // Return the generated window data.
    }

    // ------------------------------------ //
    // Window Definition Methods
    // Reference: https://en.wikipedia.org/wiki/Window_function
    // https://web.archive.org/web/20050113013738id_/http://congres.cran.uhp-nancy.fr/ICASSP_2001/MAIN/papers/pap45.pdf
    // https://www.ni.com/docs/en-US/bundle/labwindows-cvi/page/advancedanalysisconcepts/lvac_low_sidelobe.html?srsltid=AfmBOoq24bE811jsNCA5Frywall7E4fABxA6kj3FgSxqYY_808W37dA1

    // ------------------------------------ //

    inline vector<T> Hanning(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
            w[n] = 0.5 * (1 - cos(2 * M_PI * n / (N - 1)));
        return w;
    }


    inline vector<T> Hamming(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
            w[n] = 0.5383553946707251 - 0.4616446053292749 * cos(2 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> BlackmanHarris(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
            w[n] = 0.35875 - 0.48829 * cos(2 * M_PI * n / (N - 1)) + 0.14128 * cos(4 * M_PI * n / (N - 1)) - 0.01168 * cos(6 * M_PI * n / (N - 1));
        return w;
    }

    inline vector<T> ExactBlackman(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
            w[n]= .4243800934609435 - 0.4973406350967378 * cos(2 * M_PI * n / (N - 1)) + 0.7827927144231873 * cos(4 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> Blackman(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
            w[n] = 0.42 - 0.5 * cos(2 * M_PI * n / (N - 1)) + 0.08 * cos(4 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> FlatTop(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
            w[n] = 0.21557895 - 0.41663158 * cos(2 * M_PI * n / (N - 1)) + 0.277263158 * cos(4 * M_PI * n / (N - 1)) - 0.083578947 * cos(6 * M_PI * n / (N - 1)) + 0.006947368 * cos(8 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> FourTermBHarris(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
        {
            w[n] = 0.3635819267707608 - 0.4891774371450171 * cos(2 * M_PI * n / (N - 1)) + 0.1365995139786921 * cos(4 * M_PI * n / (N - 1)) - 0.01064112210553003 * cos(6 * M_PI * n / (N - 1));
        }
        return w;
    }


    inline vector<T> SevenTermBHarris(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
        {
            w[n] = 0.27105140069342 - 0.43329793923448 * cos(2 * M_PI * n / (N - 1)) + 0.21812299954311 * cos(4 * M_PI * n / (N - 1)) - 0.06592544638803 * cos(6 * M_PI * n / (N - 1)) + 0.01081174209837 * cos(8 * M_PI * n / (N - 1)) - 0.00077658482522 * cos(10 * M_PI * n / (N - 1)) + 0.00001388721735 * cos(12 * M_PI * n / (N - 1));
        }
        return w;
    }


    inline vector<T> LowSideLobe(const size_t N)
    {
        vector<T> w(N, T(0));
        for (size_t n = 0; n < N; ++n)
        {
            w[n] = 0.471492057 - 0.17553428 * cos(2 * M_PI * n / (N - 1)) + 0.028497078 * cos(4 * M_PI * n / (N - 1)) - 0.001261367 * cos(6 * M_PI * n / (N - 1));
        }
        return w;
    }
    inline vector<T> Tukey(
        const size_t N,                    // Length of the window.
        T alpha = 0.5)       // Proportion of the window that is tapered using a cosine taper.
    {
        vector<T> w(N, T(0));
      if (N<=0)                         // Do we have a length?
        return w;                       // If N is zero or negative, return empty vector.
      if (alpha<0)                 // Is alpha less than zero?
        alpha=0;                        // Set alpha to zero.
      if (alpha>1)                 // Is alpha greater than one?
        alpha=1;                        // Set alpha to one.
      {
        w[0]=T(1);                      // If N is 1, return a vector with one element set to 1.
        return w;                       // Return the window.
      }
      size_t taperLength=static_cast<size_t>(N*alpha/2.0);
      for (size_t i=0;i<N;++i)          // Loop through the samples.
      {
        w[i]=0.5*(1+std::cos(M_PI)*(-1+(2.0*i)/(alpha+N))); // Calculate the taper.
        if (i<taperLength)              // Sample less than taper lengh?
          w[i]=0.5*(1+std::cos(M_PI*(-1+(2.0*i)/(alpha*N))));
        else if (i>=N-taperLength)     // Sample greater than or equal to N-taperlength?
          w[i]=0.5*(1+std::cos(M_PI*(1-(2.0*(N-1-i))/(alpha*N))));
        else
          w[i]=T(1);              // Otherwise, set the value to 1.  
      }
        return w;
    }
    inline vector<T> Bartlett(const size_t N)
    {                                   // ----------- Bartlett ------------------ //
      std::vector<T> w(N,T(0));         // Create a vector of size N initialized to zero.
      if (N<=0)                         // Do we have a length?
        return w;                       // If N is zero or negative, return empty vector.
      else if (N==1)                    // Just one tap?
      {
        w[0]=T(1);                      // If N is 1, return a vector with one element set to 1.
        return w;                       // Return the Bartlett window.
      }
      size_t M=N-1;                            //
      for (size_t i=0;i<M;++i)          // Loop throught the samples.
      {                                 // and Bartlett window.
        if (i<=M/2)                     // Less or equal to half?
          w[i]=2.0*i/M;                 // Set the value to this.
        else
          w[i]=2.0-(2.0*i/M);           
      }                                 // Set the window values.
      return w;                         // Return the Bartlett window.
    }                                   // ----------- Bartlett ------------------ //
    inline vector<T> Gaussian(
        const size_t N,                    // The length of the window.
        T sigma=0.2f)             // Width and sidelobe characterstics of Gaussian window.
    {                                   // ----------- Gaussian ------------------ //
      vector<T> w(N,T(0));              // The vector for our window.
      if (N<=0)                         // Do we have a length?
        return w;                       // If N is zero or negative, return empty vector.
      else if (N==1)                    // Just one tap?
      {
        w[0]=T(1);                      // If N is 1, return a vector with one element set to 1.
        return w;                       // Return the window.
      }                           
      T M=(N-T(1))/2.0f;                // Calculate the center of the window.
      for (size_t n=0;n<N;++n)             // Loop through the length of the window.
      {                                 // and generate the window.
        T ex=-0.5*std::pow((n-M)/sigma,2.0);// Euler's exponent.
        w[n]=std::exp(ex);              // The window val at this sample.
      }                                 // Done calculating the window.
      return w;                         // Return the Gaussian window.
    }                                   // ----------- Gaussian ------------------ //
    inline vector<T> Kaiser(
      size_t N,                         // The length of the window.
      T sigmaK=4.0f)                    // The sigma parameter for the Kaiser window.
    {                                   // ------------- Kaiser ------------------ //
      vector<T> w(N,T(0));              // The vector for our window.
      if (N<=0)                         // Do we have a length?
        return w;                       // If N is zero or negative, return empty vector.
      else if (N==1)                    // Just one tap?
      {
        w[0]=T(1);                      // If N is 1, return a vector with one element set to 1.
        return w;                       // Return the window.
      }
      T beta=M_PI*sigma;                // Calculate Kaiser's beta parameter.
      T i0_beta=besselI0(beta);         // The I0th modified bessel function term.
      T alpha=(N-1)/2.0f;               // The alpha term of the window.
      for (size_t n=0;n<N;++n)          // For each point of the window...
      {
        T term=(n-alpha)/alpha;         // Get the square rooted term.
        w[n]=besselI0(beta*std::sqrt(1-term*term))/i0_beta;
      }
      return w;                         // Return the Kaiser window.
    }                                   // ------------- Kaiser ------------------ //
    inline vector<T> Rectangular(const size_t N)
    {
        vector<T> w(N, T(1)); // Rectangular window is all ones.
      if (N<=0)                         // Do we have a length?
        return w;                       // If N is zero or negative, return empty vector.
      else if (N==1)                    // Just one tap?
      {
        w[0]=T(1);                      // If N is 1, return a vector with one element set to 1.
        return w;                       // Return the window.
      }                     
        return w;
    }

   
private:
    size_t windowsize;
    T alpha{0.5}; // Default alpha for Tukey window.
    T sigma{0.2}; // Gaussian lobe width.
    WindowType window;
    vector<T> data;
    // Function to calculate the modified Bessel function of the first kind.
    T besselI0(T x)
    {                                   // ----------- besselI0
      T sum=T(1);                       // Where to store the sum.
      T term=T(1);                      // The term in the series.
      size_t k=1;                       //  Loop counter.
      while (term>1e-15f)               // While our term is greater than almost zero.
      {
        term*=(x/2.0f)*(x/2.0f)/(k*k);  // Calculate the next term in the series.
        sum+=term;                      // Add the term to the sum.
        k++;                            // Increment the loop counter.
      }                                 //
      return sum;                       // Return calculated Bessel function.   
    }
};

}

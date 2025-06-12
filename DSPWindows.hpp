#pragma once

#include <vector>
#include <complex>
#include <cmath>
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
        Rectangular
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
    // Accessors
   inline  const Window<T> GetWindow (void) const { return window; }
   const inline int GetWindowsize (void) const { return windowsize; }
   inline vector<T> GetDefaultWindow (void) { return Rectangular(windowsize);}  
   inline void SetWindowsize (const int wsiz) {windowsize=wsiz;}
   inline WindowType GetWindowType (void) const { return window; }   
    // ------------------------------------ //
    // Constructors and Destructors
    // ------------------------------------ //
    Window(const int N)
        : windowsize(N), window(WindowType::Rectangular) {}
    Window(const WindowType &w, const int N) 
    {
    SetWindowType(w,N);
    }

    inline void SetWindowType(const WindowType &w, const int N)
    {
    window=w;                            // Store the type
    windowsize=N;                        // This long.
    GenerateWindow(w, N);               // Generate the window.
    }
    inline vector<T> GenerateWindow(const WindowType& w, const int N)
    {
        switch (w)                          // Set windows according to window type.
        {
            case WindowType::Hanning:         data=Hanning(N);
            case WindowType::Hamming:         data=Hamming(N);
            case WindowType::BlackmanHarris:  data=BlackmanHarris(N);
            case WindowType::ExactBlackman:   data=ExactBlackman(N);
            case WindowType::Blackman:        data=Blackman(N);
            case WindowType::FlatTop:         data=FlatTop(N);
            case WindowType::FourTermBHarris: data=FourTermBHarris(N);
            case WindowType::SevenTermBHarris:data=SevenTermBHarris(N);
            case WindowType::LowSideLobe:     data=LowSideLobe(N);
            case WindowType::Rectangular:     data=Rectangular(N);
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

    inline vector<T> Hanning(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
            w[n] = 0.5 * (1 - cos(2 * M_PI * n / (N - 1)));
        return w;
    }


    inline vector<T> Hamming(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
            w[n] = 0.5383553946707251 - 0.4616446053292749 * cos(2 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> BlackmanHarris(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
            w[n] = 0.35875 - 0.48829 * cos(2 * M_PI * n / (N - 1)) + 0.14128 * cos(4 * M_PI * n / (N - 1)) - 0.01168 * cos(6 * M_PI * n / (N - 1));
        return w;
    }

    inline vector<T> ExactBlackman(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
            w[n]= .4243800934609435 - 0.4973406350967378 * cos(2 * M_PI * n / (N - 1)) + 0.7827927144231873 * cos(4 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> Blackman(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
            w[n] = 0.42 - 0.5 * cos(2 * M_PI * n / (N - 1)) + 0.08 * cos(4 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> FlatTop(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
            w[n] = 0.21557895 - 0.41663158 * cos(2 * M_PI * n / (N - 1)) + 0.277263158 * cos(4 * M_PI * n / (N - 1)) - 0.083578947 * cos(6 * M_PI * n / (N - 1)) + 0.006947368 * cos(8 * M_PI * n / (N - 1));
        return w;
    }


    inline vector<T> FourTermBHarris(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
        {
            w[n] = 0.3635819267707608 - 0.4891774371450171 * cos(2 * M_PI * n / (N - 1)) + 0.1365995139786921 * cos(4 * M_PI * n / (N - 1)) - 0.01064112210553003 * cos(6 * M_PI * n / (N - 1));
        }
        return w;
    }


    inline vector<T> SevenTermBHarris(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
        {
            w[n] = 0.27105140069342 - 0.43329793923448 * cos(2 * M_PI * n / (N - 1)) + 0.21812299954311 * cos(4 * M_PI * n / (N - 1)) - 0.06592544638803 * cos(6 * M_PI * n / (N - 1)) + 0.01081174209837 * cos(8 * M_PI * n / (N - 1)) - 0.00077658482522 * cos(10 * M_PI * n / (N - 1)) + 0.00001388721735 * cos(12 * M_PI * n / (N - 1));
        }
        return w;
    }


    vector<T> LowSideLobe(const int N)
    {
        vector<T> w(N, T(0));
        for (int n = 0; n < N; ++n)
        {
            w[n] = 0.471492057 - 0.17553428 * cos(2 * M_PI * n / (N - 1)) + 0.028497078 * cos(4 * M_PI * n / (N - 1)) - 0.001261367 * cos(6 * M_PI * n / (N - 1));
        }
        return w;
    }


    vector<T> Rectangular(const int N)
    {
        vector<T> w(N, T(1));
        return w;
} 

   
private:
    int windowsize;
    WindowType window;
    vector<T> data;
};

}

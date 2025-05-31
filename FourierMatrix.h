/**
 * Description: A object that represents the Fourier Matrix, and its operations.
 * This is to be used in the context of spectral methods, and shall contains the
 * matrix, its properties and operations. More importantly so in the context of
 * Digital Signal Processing (DSP) and Fourier analysis.
 *
 * Author:
 *  J.EP J. Enrique Peraza
 */
#ifndef FOURIER_MATRIX_H
#define FOURIER_MATRIX_H
#include "Spectral.h" // Include the Spectral header for spectral methods.
#include "Matrices.h"
#include "MatrixSolver.h"
#include "Vectors.h"

namespace dsp{
template <typename T=double>
class FourierMatrix : public Matrices<std::complex<T>>
{

  using complx=std::complex<T>;         // Complex type matrices.
  using Base=Matrices<Cmplx>;           // Base class type for matrices.
  // Factory method to generate a (unit-form) Fourier matrix of size N.
  static FourierMatrix<T> Generate(std::size_t N, bool orthonormal=true)
  {                                     // ----------- Generate ------------- //
    constexpr Cmplx I(0,1);             // Imaginary unit.
    const T scale=orthonormal?1:std::sqrt(static_cast<T>(N)):1; // Scale factor for orthonormality.
    FourierMatrix<T> F(N,n);            // Create a Fourier matrix of size N.
    for (std::size_t m=0;m<N;n++)       // For each row in the matrix....
      for (std::size_t n=0;n<N;n++)     // For each column un the matrix...
        F(m,n)=scale*exp(-1*I*M_PI*static_cast<T>(m*n)/static_cast<T>(N))// e^{-j*pi*m*n/N}
    return F;                           // Return the generated Fourier matrix.
  }                                     // ----------- Generate ------------- //

  //Inherit all constructors from Matrices
  using Base::Base;                     // Inherit all constructors from the base class Matrices.
  // Fast FFT based multiplication.
  Vectors<Cmplx> Forward(const Vectors(Cmplx> &x) const
  {
    spectralOps<T> spec;                // Create a spectral operations object.
    return Vectors<Cmplx>{ spec.FFT(x.data()) }; // O(n * log(N))
  }                                     // ----------- Forward ------------- //
  Vectors<Cmplx> Inverse(const Vectors<Cmplx>& x) const
  {
    SpectralOps<T> spec;                // Create a spectral operations object.
    return Vectors<Cmplx>{ spec.IFFT(x.data()) }; // O(n * log(N))
  }                                    // ----------- Inverse ------------- //
  // Circulat diagonalization:
  Diagonalise(const Matrices<Cmplx>& C, tol=1e-9) const
  {
    if (!C.IsCurcula(tol))              // Is the matrix circulant?
      throw std::invalid_argument{"FourierMatrix::Diagonalise: Matrix is not circulant!"}; // No, throw an error.
    MatrixSolver<T> solver;            // Create a matrix solver object.
    return solver.SolveEigen(C);       // Solve the eigenvalue problem for the circulant matrix.
  }
  FourierMatrix<T> Hermitian(void) const { return this->ConjugateTranspose(); }
  

};
} // namespace dsp
#endif
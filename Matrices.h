/**
 * A class to perform matrix operations. 
 * 
 * Author:
 *  J.EP J. Enrique Peraza
 */
#ifndef MATRICES_H
#define MATRICES_H
#include <vector>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <complex>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <utility>


template <typename T>
class Matrices{
  
public:
  // Constructor and destructors
  Matrices(void)=default;               // Default constructor for an empty matrix.
  Matrices(size_t rows, size_t cols)    // Overloaded constructor
    :rows{rows},cols{cols}, mat(rows, std::vector<T>(cols, T{})){}
  Matrices(std::initializer_list<std::vector<T>> init) // Initializer list constructor
  {
    rows=init.size();                   // Set the number of rows from the initializer list size.
    cols=rows>0?init.begin()->size():0;  // Set the number of columns from the first row size, or 0 if empty.
    mat=std::vector<std::vector<T>>(init);// Initialize the matrix with the init list.
    // validate uniform column count:
    for (const auto& row:mat)           // Loop through each row in the matrix.
      if (row.size()!=cols)             // Check if the row size matches the expected column size.
        throw std::invalid_argument{"Matrices::Initializer list rows must have the same number of columns!"};
  }
  Matrices(const Matrices&)=default;    // Copy constructor.
  Matrices(Matrices&&)=default;         // Move constructor.
  virtual ~Matrices()=default;          // Default destructor.

  // Element accessors:
  // Access element at (i, j) with bounds checking (mutable).
  T& at(size_t i, size_t j) 
  {
    if (i >= mat.size() || j >= mat[i].size())
      throw std::out_of_range{"Matrices::at(): Index out of range!"};
    return mat[i][j];
  }
  // Access element at (i, j) with bounds checking (const).
  const T& at(size_t i, size_t j) const 
  {
    if (i >= mat.size() || j >= mat[i].size())
      throw std::out_of_range{"Matrices::at(): Index out of range!"};
    return mat[i][j];
  }
  // Fast element acces with no bounds checking (mutable).
  T& operator()(size_t i,size_t j) { return mat[i][j]; }
  // Fast element access with no bounds checking (const).
  const T& operator()(size_t i,size_t j) const { return mat[i][j]; }
  // rows and cols dimensional accessors:
  size_t Rows() const noexcept { return rows; } // Returns the number of rows in the matrix.
  size_t Cols() const noexcept { return cols; } // Returns the number of columns in the matrix.

  // Utility methods:
  // Fill the matrix with a value v.
  void fill(const T& v)
  {
    for (auto& row : mat)               // Loop through each row in the matrix.
      std::fill(row.begin(), row.end(), v);// Fill each row with the value v.
  }
  // Clear the matrix, removing all elements.
  void clear() noexcept 
  {
    mat.clear();                        // Clear the matrix data.
    rows = 0;                           // Reset the number of rows to 0.
    cols = 0;                           // Reset the number of columns to 0.
  }
  // Resize the matrix to new dimensions (newRows, newCols).
  void resize(size_t newRows, size_t newCols, const T& value = T{})
  {
    mat.resize(newRows, std::vector<T>(newCols, value)); // Resize the matrix to new dimensions.
    rows = newRows;                     // Update the number of rows.
    cols = newCols;                     // Update the number of columns.
  }
  // Get the underlying matrix data.
  const std::vector<std::vector<T>>& data() const noexcept { return mat; } // Returns a const reference to the matrix data.
  std::vector<std::vector<T>>& data() noexcept { return mat; } // Returns a reference to the matrix data.
  // Check if the matrix is empty.
  bool empty() const noexcept { return mat.empty(); } // Returns true if the matrix is empty.
  // Get the number of elements in the matrix.
  size_t size() const noexcept { return rows * cols; } // Returns the total number of elements in the matrix.
  // Get the number of elements in a specific row.
  size_t row_size(size_t i) const 
  {
    if (i >= rows)
      throw std::out_of_range{"Matrices::row_size(): Row index out of range!"};
    return mat[i].size();              // Returns the number of columns in the specified row.
  }
  // Get the number of elements in a specific column.
  size_t col_size(size_t j) const 
  {
    if (j >= cols)
      throw std::out_of_range{"Matrices::col_size(): Column index out of range!"};
    return std::count_if(mat.begin(), mat.end(),
      [j](const std::vector<T>& row) { return j < row.size(); }); // Count the number of rows that have at least j elements.
  }
  // Compute the Frobenius norm of the matrix. The Frobenius norm is defined as
  // the square root of the sum of the squares of all elements in the matrix.
  // It is a measure of the "size" of the matrix.
  T FrobeniusNorm() const
  {                                     // ----------- FrobeniusNorm ------------ //
    if (mat.empty() || mat[0].empty())  // Check if the matrix is empty.
      return T{};                       // Return zero if the matrix is empty.
    T sum{};                            // Initialize the sum to zero.
    size_t N=rows;                      // Get the number of rows in the matrix.
    for (size_t i=0;i<N;i++)            // Loop through each row in the matrix.
      for (const auto& elem : mat[i])   // Loop through each element in the row.
        sum+=elem*elem;                 // Add the square of the element to the sum.
    return std::sqrt(sum);              // Return the square root of the sum as the Frobenius norm.
  }                                     // ----------- FrobeniusNorm ------------ //
  // Get the average power of all signals. The row or col vectors
  // in the matrix representing each a signal with the elements representing the samples
  // of the signals.
  T FrobeniusAveragePower() const
  {
    if (mat.empty() || mat[0].empty())   // Do we have a matrix of signals to operate on?
      return T{};                        // No, can't do much without data. Exit.
    T sum{};                             // Initialize sum to zero.
    size_t N=rows;                       // Get the number of rows in the matrix.
    for (size_t i=0;i<N;i++)             // Loop through each row in the matrix.
    {
      T row_energy{};                    // The energy per row.
      for (const auto& elem:mat[i])      // For every element in the matrix.
        row_energy+=elem*elem;           // Add the square of the element to the result.
      sum+=row_energy/static_cast<T>(cols);// Aggregate the energy off the signals.
    }                                    // Done aggregating sum of energy in signal arremsbly
    return sum/static_cast<T>(rows*cols);// Divide energy by all elements of the matrix.
  }                                      // ------- FrobeniusAveragePower ----------- //  
  // Copy and move assignment operators:
  Matrices& operator=(const Matrices& rhs) = default; // Copy assignment operator.
  Matrices& operator=(Matrices&& rhs) noexcept = default; // Move assignment operator.
  // Arithmetic operations:
  Matrices operator+(const Matrices& rhs) const
  {
    if (rows!=rhs.Rows() || cols!=rhs.Cols())
      throw std::invalid_argument{"Matrices::operator+: Size mismatch: Matrices must have the same dimensions for addition!"};
    Matrices result(rows, cols);        // Create a new matrix to hold the result.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
        result.mat[i][j]=mat[i][j]+rhs.mat[i][j]; // Add the corresponding elements.
    return result;                      // Return the resulting matrix.
  }
  Matrices operator-(const Matrices& rhs) const
  {
    if (rows!=rhs.Rows() || cols!=rhs.Cols())
      throw std::invalid_argument{"Matrices::operator-: Size mismatch: Matrices must have the same dimensions for addition!"};
    Matrices result(rows, cols);        // Create a new matrix to hold the result.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
        result.mat[i][j]=mat[i][j]-rhs.mat[i][j]; // Subs the corresponding elements.
    return result;                      // Return the resulting matrix.
  }
  // Matrix multiplication.
  Matrices operator*(const Matrices& rhs) const
  {
    if (cols!=rhs.Rows())                 // Right inner dimensions for mult?
      throw std::invalid_argument{"Matrices::operator*: Size mismatch: Inner dimensions must match for multiplication!"};
    Matrices result(rows, rhs.Cols());    // Create a new matrix to hold the result.
    for (size_t i=0;i<rows;i++)         // For each row of left matrix...
      for (size_t j=0;j<rhs.Cols();j++) // For each column of right matrix...
      {                                 // Loop through elements to compute product.
        T sum{};                        // Initialize the sum for the dot product.
        for (size_t k=0;k<cols;k++)     // Loop through the inner dimension.
          sum+=mat[i][k]*rhs.mat[k][j]; // Compute the dot product.
        result.mat[i][j]=sum;           // Store the result in the new matrix.
      }                                 // Done computing the product.
    return result;                      // Return the resulting matrix.
  }
  // Scalar multiplication.
  Matrices operator*(const T& c) const
  {
    Matrices result{rows,cols};         // Create a new matrix to hold the result.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix.
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
        result.mat[i][j]=mat[i][j]*c;   // Scale each element by the scalar value.
    return result;                      // Return the resulting matrix.
  }
  // Matrix * vector multiplication.
  Matrices operator*(const std::vector<T>& vec) const
  {
    if (cols!=vec.size())               // Check if the vector size matches the number of columns.
      throw std::invalid_argument{"Matrices::operator*: Size mismatch: Vector size must match the number of columns!"};
    Matrices result(rows, 1);           // Create a new matrix to hold the result (column vector).
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
    {                                   // Compute the dot product with the vector.
      T sum{};                          // Initialize the sum for the dot product.
      for (size_t j=0;j<cols;j++)       // Loop through the columns.
        sum+=mat[i][j]*vec[j];          // Compute the dot product.
      result.mat[i][0]=sum;             // Store the result in the new matrix.
    }                                   // Done computing the product.
    return result;                      // Return the resulting matrix.
  }
  // Scalar division.
  Matrices operator/(const T& scalar) const
  {
    if (scalar==T{})                   // Check for division by zero.
    throw std::invalid_argument{"Matrices::operator/: Division by zero!"};
    Matrices result{rows,cols};        // Create a new matrix to hold the result.
    for (size_t i=0;i<rows;i++)        // For each row in the matrix.
    for (size_t j=0;j<cols;j++)        // For each column in the matrix...
        result.mat[i][j]=mat[i][j]/scalar; // Scale each element by the scalar value.
    return result;                     // Return the resulting matrix.
  }

  // Transpose:
  Matrices Transpose() const
  {
    Matrices result{cols,rows};         // Create a new matrix with transposed dimensions.
    for (size_t i=0;i<rows;i++)         // For each row in the original matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the original matrix...
        result.mat[j][i]=mat[i][j];     // Assign the transposed element.
    return result;                      // Return the transposed matrix.
  }
  /// Matrix properties assertions:
  // IsDiagonal: Check if the matrix is exactly diagonal (off-diagonal entries are 0).
  // To help us guard against trivial cases, (e.g. to skip heavy eigen-solver 
  // computations) such as the matrix already being diagonal.
  bool IsDiagonal(double tol=1e-9) const
  {                                     // ----------- IsDiagonal ------------ //
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
        if (i!=j&&std::abs(mat[i][j]>tol))// Are the off diagonals zero?
          return false;                 // No, the matrix is not diagonal.
    return true;                        // Yes, the matrix is diagonal.
  }                                     // ----------- IsDiagonal ------------ //

  // IsSymmetric: Check if the matrix is symmetric (A == A^T for real matrices)
  // or Hermitian (A == A^H for complex matrices). Permits a dispatcher pick
  // Jacobi or tridiagonalization+QL algorithm for eigenvalue computations
  // when appropiate.
  bool IsSymmetric(double tol=e-9) const
  {                                     // ----------- IsSymmetric ------------ //
    if (rows!=cols)                     // Is the matrix is square?
      return false;                     // No, it cannot be symmetric.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
      {                                 // Test for symmetry.
        auto a=mat[i][j];               // Matrix element at (i,j).
        auto b=std::is_floating_point_v<T>?mat[j][i]:std::conj(mat[j][i]); // Conjugate for complex types.
        if (std::abs(a-b)>tol)          // Are the elements equal within tolerance?
          return false;                 // No, the matrix is not symmetric.
      }                                 // Done checking all elements.
    return true;                        // Else, the matrix is symmetric.
  }                                     // ----------- IsSymmetric ------------ //
  // IsNormal: Check if the matrix is normal (A*A^H == A^H*A). Tells us if the 
  // matrix commutes with its adjoint, helpful for fast analytic solvers.
  bool IsNormal(double tol=1e-9) const
  {                                     // ----------- IsNormal ------------ //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so its not normal either.
    Matrices<T> AH=this->conjugateTranspose(); // Get the conjugate transpose of the matrix.
    Matrices<T> AAH=(*this)*AH;         // Compute A * A^H.
    Matrices<T> AHA=AH*(*this);         // Compute A^H * A.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
        if (std::abs(AAH(i,j)-AHA(i,j)>tol)) // Are the elements equal within tolerance?
          return false;                 // No, the matrix is not normal.
    return true;                        // Else, the matrix is normal.
  }                                     // ----------- IsNormal ------------ //

  // IsUnitary: Check if the matrix is unitary (A*A^H == I).
  bool IsUnitary(double tol=1e-9) const
  {                                     // ----------- IsUnitary ------------ //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so its not unitary either.
    Matrices<T> AH=this->conjugateTranspose(); // Get the conjugate transpose of the matrix.
    Matrices<T> I=(*this)*AH;           // Compute A * A^H.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
    {                                   //    and..
      for (size_t j=0;j<cols;j++)       // ...for each column in the matrix...
      {                                 // Test for unitary matrix.
        if (i==j)                       // If we are on the diagonal...
          if (std::abs(I(i,j)-T{1})>tol) // Is it equal to 1 within tolerance?
            return false;               // No, so its not unitary.
        else if (std::abs(I(i,j))>tol)  // Else, is it zero within tolerance?
          return false;                 // No, so its not unitary.
      }                                 // Done checking all elements.
    }                                   // Done testing for unitary.
    return true;                        // Else, the matrix is unitary.
  }                                     // ----------- IsUnitary ------------ //
  // IsOrthogonal: Check if the matrix is orthogonal (A*A^T == I). Help verify
  // if the unmixing matrices in the ESPRIT algorithm are well-conditioned.
  bool IsOrthogonal(double tol=1e-9) const
  {                                     // ----------- IsOrthogonal ------------ //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so it's not orthogonal either.
    Matrices<T> MMT=this->conjugateTranspose()*(*this); // Compute A * A^T.
    Matrices<T> I(rows,cols);           // Create an identity matrix of the same size.
    I,fill(T{});                        // Fill the identity matrix with zeros.
    for (size_t i=0;i<rows;i++) I(i,i)=T{1}; // Set the diagonal elements to 1.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
        if (std::abs(MMT(i,j)-I(i,j)>tol))// Are the elemenets equal within tolerance?
          return false;                 // No, matrix is not orthogonal.
    return true;                        // Else, the matrix is orthogonal.
  }                                     // ----------- IsOrthogonal ------------ //
  // Is AntiSymmetric: Check if the matrix is anti-symmetric (A^T = -A). Helps us
  // catch purely imaginary-valued skew fields in DSP applications, such as
  // quadrature signals, Hilberts transforms, etc.
  bool IsAntiSymmetric(double tol=1e-9) const
  {                                     // ----------- IsAntiSymmetric ------------ //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so it's not anti-symmetric either.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
    {                                   //   
       if (std::abs(mat[i][j]>tol))     // Are the diagonals greater than zero?
         return false                   // Then it is not anti-symmetric.
      for (size_t j=0;j<cols;j++)       // For each column in the mtrix...
      {                                 // Test for anti-symmetry.
        auto a=mat[i][j];               // Matrix element at (i,j).
        auto b=std::is_floating_point_v<T>?-mat[i][j]:-std::conj(mat[i][j]);// Conjugate for complex types.
        if (std::abs(a-b)>tol)          // Are the elements equal within tolerane?
          return false;                 // No, the matrix is not anti-symmetric.
      }                                 // Done checking all elements.
    }                                   // Done testing for anti-symmetry.
    return true;                        // Else, the matrix is anti-symmetric.          
  }                                     // ----------- IsAntiSymmetric ------------ //

 // IsHermitian: Check if the matrix is Hermitian (A == A^H for complex matrices).
  bool IsHermitian(double tol=1e-9) const
  {                                     // ----------- IsHermitian ------------ //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so it's not Hermitian either.
    for (size_t i=0;i<rows;i++)         // For each row in the matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the matrix...
      {                                 // Test for Hermitian property.
        auto a=mat[i][j];               // Matrix element at (i,j).
        auto b=std::conj(mat[j][i]);    // Conjugate for complex types.
        if (std::abs(a-b)>tol)          // Are the elements equal within tolerance?
          return false;                 // No, the matrix is not Hermitian.
      }                                 // Done checking all elements.
    return true;                        // Else, the matrix is Hermitian.
  }                                     // ----------- IsHermitian ------------ //

  // Is Circular: Check if the matrix is circulant matrix. A circular matrix
  // is a special type of Toeplitz matrix where each row is a cyclic shift of the previous row.
  // It is useful in signal processing applications, such as cyclic convolution 
  // and filtering. In other words, when multiplied by the Fourier matrix, the
  // circulant matrix gets diagonalized. Where the diagonals are related to the
  // eigenvalues of the circulant matrix. Such that for any Circulant Matrix C
  // there exists a Fourier Matrix F, such that F^-1*C*F is a diagonal matrix,
  // where the resulting diagonal elements are the eigenvalues of C.
  // As a sidebar, the columns of the Fourier matrix are the eigenvalues of the
  // circulant matrix.
  bool IsCircular(double tol=1e-9) const
  {                                     // ----------- IsCircular ----------- //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so it's not circulant either.
    // We store the first row of the matrix to compare with the rest.
    const Sstd::vector<T> &firstrow=mat[0];// Store the first row of the matrix.
    for (size t i=1;i<rows;i++)         // For each subsequent row...
    {                                   // Start testing for circulance property.
      bool foundmatch=false;            // Flag to check if we found a match.
      for (size_t shift=0;shift<cols;shift++)// For each possible shift...
      {                                 // 
        std::vector<T> shiftedrow(cols);// Build shifted version of mat[i].
        for (size_t j=0;j<cols;j++)     // For each column in the row...
          shiftedrow[j]=mat[i][(j+shift)%cols];// Shift elements cyclically.
        if constexpr (std::is_floating_pint_v<T>)// Is it a floating point vector?
        {                               // Yes so we should factor in tolerance.
          bool allclose=true;           // Flag to check if all elements are close.
          for (size_t j=0;j<cols;j++)   // For each column in the row....
          {                             // Check if we are within tolerance...
            if (std::abs(shiftedrow[j]-firstrow[j])>tol)// Are the elements close?
            {                           // No..
              allclose=false;           //  ..set the flag to false.
              break;                    // Break out of the test loop.
            }                           // Done checking all elements.
          }                             // Done checking all elements in the row.
          if (allclose)                 // if all element are close...
          {                             // We found a match.
            foundmatch=true;            // Set the flag to true.
            break;                      // Break out of the shift loop.
          }                             // Done checking if all shift elems are close.  
        }                               // Done checking if floaing point.
        else                            // Else, not flaoting point.
        {                               // So we can check for exact equality.
          if (shiftedrow==firstrow)     // Are the shifted and first row equal?
          {                             // Yes, found a match.
            foundmatch=true;            // Set the flag to true.
            break;                      // Break out of the shift loop.
          }                             // Done checking if shifted row is equal to first row.
        }                               // Done checking if not floating point.
      }                                 // Done checking all shifts.
      if (!foundmatch)                  // If we did not find a match...
        return false;                   // Then the matrix is not circulant.
    }                                   // Done checking all rows.
    return true;                        // If we got here the matrix is circulant.
  }                                     // ----------- IsCircular ----------- // 
  // IsToeplitz: Check if a matrix is a Toeplitz matrix. A toeplitz matrix is constant
  // along its diagonals. So we check if A[i][j]=A[i-1][j-1] for all i,j.
  // Topelitz matrix in DSP represent Linear Time Invariant (LTI) systems,
  // where the system response is constant over time. This property is useful
  // in signal processing applications, such as filtering and convolution.
  bool IsToeplitz(double tol=1e-9) const
  {                                     // ----------- IsToeplitz ------------ //
    if (rows==0||cols==0)               // Empty matrix?
      return false;                     // Can't do much.
    for (size_t i=1;i<rows;i++)         // For the first row...
    {                                   //  We start checking the diagonals.
      for (size_t j=1;j<cols;j++)       // For the first column...
      {                                 // 
        if constexpr (std::is_floating_point_v<T>)// Is it a floating point matrix?
        {                               // Yes, so we should factor in tolerance.
          if (std::abs(mat[i][j]-mat[i-1][j-1])>tol) // Are the elements equal within tolerance?
            return false;               // No, so it's not Toeplitz.
        }                               // Done checking if floating point.
        else                            // Else, not floating point.
        {                               // So we can check for exact equality.
          if (mat[i][j]!=mat[i-1][j-1]) // Are the elements equal?
            return false;               // No, so it's not Toeplitz.
        }                               // Done checking if not floating point.
      }                                 // Done checking all columns in the row.
    }                                   // Done checking all rows.
    return true;                        // If we got here, the matrix is Toeplitz.
  }                                     // ----------- IsToeplitz ------------ //
  // IsPositiveDefinite: Check if the matrix is positive definite.
  // A positive definite matrix is a symmetric matrix with all positive eigenvalues.
  // It attempts an in-place Cholesky. If any pivot is less or equal to tol, then
  // it is not positive definite. It is true if A is strictly (Hermitian) positive definite.
  // 
  bool IsPositiveDefinite(double tol=1e-9) const
  {                                     // ----------- IsPositiveDefinite ------------ //
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so it's not positive definite either.
    Matrices<T> C=*this;                // Create a copy so we do not overwite original.
    for (size_t k=0;k<rows;k++)         // For each row in the matrix...
    {                                   // Begin computations...
    // Compute the sum{m=0...k-1} |C(k,m)|^2
      long double sum=0.0L;             // Initialize the sum to zero.
      for (size_t m=0;m<rows;m++)       // For each row in the matrix.
      {                                 // Start testing matrix.
        if constexpr (std::is_floating_pint_v<T>)// Are elems floating point?
        {                               // Yes so get the norm.
          sum+=std::norm(C(k,m));       // Compute the norm of the element.
        }                               // Done checking if floating point.
        else                            // Else, not floating point.
        {                               // Se we can treat norm as squared.
          sum+=static_cast<long double>(C(k,m))*
            static_cast<long double>(C(k,m)); // Compute the squared element.
        }                               // Done checking if not floating point.
      }                                 // Done computing the sum.
      long double diag=static_cast<long double>(C(k,k))-sum; // Compute the diagonal element.
      if (diag<=tol)                    // Is pivot less than or equal to tol?
        return false;                   // Yes, so the matrix is not positive definite.
      C(k,k)=static_cast<T>(std::sqrt(diag));// Set the diagonal to sqrt(diag).
      // ------------------------------ //
      // Now compute Cholesky factorization where:
      // L(i,k)=[C(i,k)-sum{m<k} L(i,m)*L(k,m)]/L(k,k)
      // ------------------------------ //
      for (size_t i=k+1;i<rows;i++)     // For each row below the pivot...
      {                                 // Start computing the Cholesky factorization.
        long double sum2{0.0L};         // Initialize the sum to zero.
        for (size_t m=0;m<k;m++)        // For each col in the row...
        {                               // Start computing the sum of products.
          if constexpr (std::is_floating_point_v<T>)// Are elems floating point?
          {                             // Yes,
            sum2+=C(i,m)*std::conj(C(k,m)); // Compute the sum of products.
          }                             // Done checking if floating point.
          else                          // Else, not floating point.
          {                             // So we can treat norm as squared.
            sum2+=static_cast<long double>(C(i,m))*
              static_cast<long double>(C(k,m)); // Compute the squared element.
          }                             // Done checking if not floating point.
        }                               // Done computing the sum of products.
        long double num=static_cast<long double>(C(i,k))-sum2; // Compute the numerator.
        long double denom=static_cast<long double>(C(k,k));// Get the denominator.
        if (std::abs(denom)<tol)        // Is the denominator less than tol (zero pivot)?
          return false;                 // Yes, we can't divide by zero.. so just exit.
        C(i,k)=static_cast<T>(num/denom); // Compute the Cholesky factorization.
      }                                 // Done computing the Cholesky factorization.
    }                                   // Done with all computations.
    return true;                        // If we got here, the matrix is positive definite.
  }                                     // ----------- IsPositiveDefinite ------------ //
  // IsPositiveSemiDefinite: Check if the matrix is positive semi-definite.
  // Similar to IsPositiveDefinite in the way it operates, but it allows zero pivots
  // within tol. True iff all eigenvalues >=0 (Hermitian semidifinite).
  bool IsSemiPositiveDefinite(double tol=1e-9) const
  {                                     // ----------- IsSemiPositiveDefinite 
    if (rows!=cols)                     // Is the matrix square?
      return false;                     // No, so it's not semi-positive definite either.
    Matrices<T> C=*this;                // Create a copy so we do not overwite original.
    for (size_t k=0;k<rows;k++)         // For each row in the matrix...
    {                                   // Begin computations...
      // Compute the sum{m=0...k-1} |C(k,m)|^2
      long double sum=0.0L;             // Initialize the sum to zero.
      for (size_t m=0;m<rows;m++)       // For each row in the matrix.
      {                                 // Start testing matrix.
        if constexpr (std::is_floating_pint_v<T>)// Are elems floating point?
        {                               // Yes so get the norm.
          sum+=std::norm(C(k,m));       // Compute the norm of the element.
        }                               // Done checking if floating point.
        else                            // Else, not floating point.
        {                               // Se we can treat norm as squared.
          sum+=static_cast<long double>(C(k,m))*
            static_cast<long double>(C(k,m)); // Compute the squared element.
        }                               // Done checking if not floating point.
      }                                 // Done computing the sum.
      long double diag=static_cast<long double>(C(k,k))-sum; // Compute the diagonal element.
      if (diag<-tol)                    // Is pivot less than negative tol?
        return false;                   // Yes, so the matrix is not positive semi-definite.
      // Allow for diag to approximate zero.
      if (diag<=tol)                    // Is pivot less than or equal to tol?
      {                                 // Yes, so we set the diagonal to zero.
        C(k,k)=T{0};                    // Set the diagonal to zero if it is within tolerance.
        continue;                       // Skip the rest of the computations for this pivot.
      }                                 // Done checking if pivot is less than or equal to tol.
       C(k,k)=static_cast<T>(std::sqrt(diag));// Set the diagonal to sqrt(diag).
      // ------------------------------ //
      // Now compute Cholesky factorization where:
      // L(i,k)=[C(i,k)-sum{m<k} L(i,m)*L(k,m)]/L(k,k)
      // ------------------------------ //
      for (size_t i=k+1;i<rows;i++)     // For each row below the pivot...
      {                                 // Start computing the Cholesky factorization.
        long double sum2{0.0L};         // Initialize the sum to zero.
        for (size_t m=0;m<k;m++)        // For each col in the row...
        {                               // Start computing the sum of products.
          if constexpr (std::is_floating_point_v<T>)// Are elems floating point?
          {                             // Yes,
            sum2+=C(i,m)*std::conj(C(k,m)); // Compute the sum of products.
          }                             // Done checking if floating point.
          else                          // Else, not floating point.
          {                             // So we
            sum2+=static_cast<long double>(C(i,m))*
              static_cast<long double>(C(k,m)); // Compute the squared element.
          }                             // Done checking if not floating point.
        }                               // Done computing the sum of products.
        long double num=static_cast<long double>(C(i,k))-sum2; // Compute the numerator.
        long double denom=static_cast<long double>(C(k,k));// Get the denominator.
        if (std::abs(denom)<tol)        // Is the denominator less than tol (zero pivot)?
          return false;                 // Yes, so we can't divide by zero...skip it.
        C(i,k)=static_cast<T>(num/denom); // Compute the Cholesky factorization.
      }                                 // Done computing the Cholesky factorization.
    }                                   // Done with all computations.
    return true;                        // If we got here, the matrix is positive semi-definite.
  }
  // Debug print.
  void Print(std::ostream& os = std::cout) const
  {                                     // ----------- Print ------------ //
    for (const auto& row:mat)           // For each row in the matrix...
    {                                   //  and...
      for (const auto& elem:row)        //   ...for each element in the row...
        os<<elem<<" ";                  // Print the element followed by a space.
      os<<std::endl;                    // Print a newline after each row.
    }                                   // Done printing the matrix.
  }                                     // ----------- Print ------------ //
private:
  // Helper function for complext transpose:
  Matrices<T> conjugateTranspose(void) const
  {                                     // --------- conjugateTranspose ----------- //
    if (mat.empty() || mat[0].empty())  // Check if the matrix is empty.
      return Matrices<T>{};             // Return an empty matrix if so.
    size_t rows=mat.size();             // Get the number of rows in the matrix.
    size_t cols=mat[0].size();          // Get the number of columns in the matrix.
    Matrices<T> C(cols,rows);           // Create a new matrix with transposed dimensions.
    for (size_t i=0;i<rows;i++)         // For each row in the original matrix...
      for (size_t j=0;j<cols;j++)       // For each column in the original matrix...
        C(j,i)=std::is_floating_pint_v<T>?mat[i][j]:std::conj(mat[i][j]); // Assign the conjugated transposed element.
    return C;                           // Return the conjugated transposed matrix.
  }                                     // --------- conjugateTranspose ----------- //
  // Data members:
    std::vector<std::vector<T>> mat;// The matrix data stored as a vector of vectors.
    size_t rows{0};                     // Number of rows in the matrix.
    size_t cols{0};                     // Number of columns in the matrix.
};
#endif

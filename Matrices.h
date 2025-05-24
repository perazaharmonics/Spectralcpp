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
  // Debug print.
  void Print(std::ostream& os = std::cout) const
  {
    for (const auto& row:mat)           // For each row in the matrix...
    {
      for (const auto& elem:row)        // For each element in the row...
        os<<elem<<" ";                  // Print the element followed by a space.
      os<<std::endl;                    // Print a newline after each row.
    }
  }
private:
    std::vector<std::vector<T>> mat;// The matrix data stored as a vector of vectors.
    size_t rows{0};                     // Number of rows in the matrix.
    size_t cols{0};                     // Number of columns in the matrix.
};
#endif
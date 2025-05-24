/**
 * A class that contains methods for vector operations. 
 * 
 * Author: 
 *   J.EP J. Enrique Peraza
 */
#ifndef VECTORS_H
#define VECTORS_H
#include <vector>
#include <complex>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <utility>

template <typename T>
class Vectors{
  using Vec = std::vector<T>;          // A vector of type T
  using Cx  = std::complex<T>;         // A complex number of type T
  using CxVec = std::vector<Cx>;       // A vector of complex numbers of type T

public:
  // Constructors and destructors:
  Vectors()=default;                    // Empty vector constructor.
  explicit Vector(size_t n):            // Zero initialized length n 
    vect(n,T{}){}                       // vector constructor.
  Vectors(std::initializer_list<T> init):// Initializer list constructor allows
    vect(init){}                        // to initialize the vector with a list of values.
  Vectors(const Vec &v):                // Copy constructor to initialize the vector with another vector.
    vect(v){}                           // This allows for easy initialization from another vector.
  Vectors(Vec &&v):                     // Move constructor to initialize the vector with another vector.
    vect(std::move(v)){}                // This allows for efficient initialization from another vector.
   virtual ~Vectors()=default;

  // Operator overriding to allow for easy manipulation of vectors:
  Vec &operator=(const Vec &v)          // Assignment operator to copy a vector.
  {                                     // ----------- operator= ------------- //
    vect=v;                             // Copy the vector.
    return vect;                        // Return the vector.
  }                                     // ----------- operator= ------------- //
  Vec &operator=(Vec &&v)               // Move assignment operator to move a vector.
  {                                     // ----------- operator= ------------- //
    vect=std::move(v);                  // Move the vector.
    return vect;                        // Return the vector.
  }                                     // ----------- operator= ------------- //
  T &operator[](size_t i)               // Subscript operator to access an element of the vector.
  {                                     // ----------- operator[] ------------- //
    if (i>=vect.size())                 // Is the index out of bounds?
    {                                   // Yes, print error message.
      std::cerr<<"[ERROR] Index out of bounds!\n";
      return vect[0];                   // Return the first element if the index is out of bounds.
    }                                   // Done checking if the index is out of bounds.
    return vect[i];                     // Return the element at the index.
  }                                     // ----------- operator[] ------------- //
  const T &operator[](size_t i) const   // Const subscript operator to access an element of the vector.
  {                                     // ----------- operator[] ------------- //
    if (i>=vect.size())                 // Is the index out of bounds?
    {                                   // Yes, print error message.
      std::cerr<<"[ERROR] Index out of bounds!\n";
      return vect[0];                   // Return the first element if the index is out of bounds.
    }                                   // Done checking if the index is out of bounds.
    return vect[i];                     // Return the element at the index.
  }                                     // ----------- operator[] ------------- //
  Vec &operator+=(const Vec &v)         // Addition assignment operator to add a vector.
  {                                     // ----------- operator+= ------------ //
    if (vect.size()!=v.size())          // Are the vectors the same length?
    {                                   // No, can't add vectors.
      std::cerr<<"[ERROR] Vectors must have the same length for addition!\n";
      return vect;                      // Return the vector if they are not the same length.
    }                                   // Done checking if the vectors are the same length.
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      vect[i]+=v[i];                    // Add the elements of the vectors.
    }                                   // Done adding the vectors.
    return vect;                        // Return the sum of the vectors.
  }                                     // ----------- operator+= ------------ //
  Vec &operator-=(const Vec &v)         // Subtraction assignment operator to subtract a vector.
  {                                     // ----------- operator-= ------------ //
    if (vect.size()!=v.size())          // Are the vectors the same length?
    {                                   // No, can't subtract vectors.
      std::cerr<<"[ERROR] Vectors must have the same length for subtraction!\n";
      return vect;                      // Return the vector if they are not the same length.
    }                                   // Done checking if the vectors are the same length.
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      vect[i]-=v[i];                    // Subtract the elements of the vectors.
    }                                   // Done subtracting the vectors.
    return vect;                        // Return the difference of the vectors.
  }                                     // ----------- operator-= ------------ //
  Vec &operator*=(const T &scalar)      // Multiplication assignment operator to scale a vector.
  {                                     // ----------- operator*= ------------ //
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      vect[i]*=scalar;                  // Scale each element by the scalar value.
    }                                   // Done scaling the vector.
    return vect;                        // Return the scaled vector.
  }                                     // ----------- operator*= ------------ //
  Vec &operator/=(const T &scalar)      // Division assignment operator to scale a vector.
  {                                     // ----------- operator/= ------------ //
    if (scalar==T(0))                   // Is the scalar zero?
    {                                   // Yes, can't divide by zero.
      std::cerr<<"[ERROR] Cannot divide by zero!\n";
      return vect;                      // Return the vector if the scalar is zero.
    }                                   // Done checking if the scalar is zero.
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      vect[i]/=scalar;                  // Scale each element by the scalar value.
    }                                   // Done scaling the vector.
    return vect;                        // Return the scaled vector.
  }                                     // ----------- operator/= ------------ //
  // ---------------------------------- //
  // Get the size of the vector.
  // ---------------------------------- //
  size_t Size() const                   // Returns the size of the vector.
  {                                     // ----------- Size ------------------ //
    return vect.size();                 // Return the size of the vector.
  }                                     // ----------- Size ------------------ //
  // ---------------------------------- //
  // Getter for the vector.
  // ---------------------------------- //
  const Vec &GetVector() const          // Returns the vector.
  {                                     // ----------- GetVector ------------ //
    return vect;                        // Return the vector.
  }                                     // ----------- GetVector ------------ //
  // ---------------------------------- //
  // Setter for the vector.
  // ---------------------------------- //
  void SetVector(const Vec &v)          // Sets the vector.
  {                                     // ----------- SetVector ------------ //
    vect=v;                             // Set the vector to the input vector.
  }                                     // ----------- SetVector ------------ //

  // Vector operations:
  // ---------------------------------- //
  // Compute the dot product of two vectors. The vector dot product tell us
  // The similarity between two vectors. The angle between the two vectors
  // can be computed from the dot product using the formula:
  // cos(theta) = (v1 . v2) / (||v1|| * ||v2||)
  // where v1 and v2 are the vectors, ||v1|| and ||v2|| are the norms of the vectors,
  // and theta is the angle between the vectors.
  // ---------------------------------- //
  T VectorDotProduct(const Vec &v1, const Vec &v2) const{
    if (v1.size()!=v2.size())           // Are they the same length?
    {                                   // No, can't compute dot product.
      std::cerr<<"[ERROR] Vectors must have the same length for dot product!\n";
      return T(0);                      // Return 0 if they are not the same length.
    }                                   // Done checking if they are the same length.
    T sum{0.0};                         // Initialize the sum.
    size_t N=v1.size();                 // Get the size of the vectors.
    for (size_t i=0;i<N;i++)            // For the length of the vectors...
    {                                   // Loop over the elements of the vectors.
      if (std::isnan(v1[i])||std::isnan(v2[i]))// Any NaN in the vectors?
      {                                 // Yes, print warning.
        std::cerr<<"[WARNING] NaN detected in vectors!\n";// Return error message.
        return T(0);                    // Return 0 if there is a NaN.
      }                                 // Done checking for NaN.
      if (std::isinf(v1[i])||std::isinf(v2[i]))// Any Inf in the vectors?
      {                                 // Yes, print warning.
        std::cerr<<"[WARNING] Inf detected in vectors!\n";
        return T(0);                    // Return 0 if there is an Inf.
      }                                 // Done checking for Inf.
        sum+=v1[i]*v2[i];               // Accumulate the product.
    }                                   // Done computing the dot product.
    return sum;                         // Return the dot product.    
  }                                     // --------- VectorDotProduct ------- //
  
  // ---------------------------------- //
  // Compute the norm of a vector. The vector norm is a measure of the length
  // of the vector. It is defined as the square root of the sum of the squares
  // of the elements of the vector. The norm is used to normalize vectors and
  // to compute the angle between vectors.
  // The norm is also known as the Euclidean norm or L2 norm. It is the equivalent
  // of the magnitute of a vector in physics. Thus it is also the power calculation.
  // ---------------------------------- //

  T VectorNorm(const Vec &v) const      // Returns norm of a vector.
  {                                     // ----------- VectorNorm ------------ //
    if (v.empty())                      // Is the vector emtpy?
    {                                   // Yes, return 0.
      std::cerr<<"[ERROR] Vector is emtpy!\n";
      return T(0);                      // Return 0 if the vector is empty.
    }                                   // Done checking if the vector is empty.
    T sum{0.0};                         // Initialize the sum.
    size_t N=v.size();                  // Get the size of the vector.
    for (size_t i=0;i<N;i++)            // For the length of the vector...
    {                                   // Loop over the elements of the vector.
      if (std::isnan(v[i]))             // Any NaN in the vector?
      {                                 // Yes print warning.
        std::cerr<<"[WARNING] NaN detected in vector!\n";
        return T(0);                    // Return 0 if there is a NaN.
      }                                 // Done checking for NaN.
      if (std::isinf(v[i]))             // Any Inf in the vector?
      {                                 // Yes, print warning.
        std::cerr<<"[WARNING] Inf detected in vector!\n";
        return T(0);                    // Return 0 if there is an Inf.
      }                                 // Done checking for Inf.
      sum+=v[i]*v[i];                   // Accumulate the square of the element.
    }                                   // Done computing the squared sum of the vector.
    return std::sqrt(sum);              // Return the square root of the sum.
  }                                     // ----------- VectorNorm ------------ //   
  
  // ---------------------------------- //
  // Transform a vector into a unit vector by dividing each element by the vector
  // norm. This is useful for normalizing vectors in spectral analysis.
  // The function returns a new vector that is the normalized version of the input vector.
  // ---------------------------------- //
  Vec NormalizeVector(const Vec &v)  const// Retuns a normalized vector.  
  {                                     // ----------- NormalizeVector --------- //
    if (v.empty())                      // Is the vector emtpy?
    {                                   // Yes, return an empty vector.
      std::cerr<<"[ERROR] Vector is emtpy!\n";
      return Vec{};                     // Return an empty vector if the vector is empty.
    }                                   // Done checking if the vector is empty.
    T norm=VectorNorm(v);               // Get the norm of the vector.
    if (norm==T(0))                     // Is the norm zero?
    {                                   // Yes, return an empty vector.
      std::cerr<<"[ERROR] Vector norm is zero!\n";
      return Vec{};                     // Return an empty vector if the norm is zero.
    }                                   // Done checking if the norm is zero.
    Vec normal(v.size());               // Create a new vector of the same size as v.
    for (size_t i=0;i<v.size();i++)     // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      normal[i]=v[i]/norm;              // Divide each element by the norm.
    }                                   // Done normalizing the vector.
    return normal;                      // Return the normalized vector.
  }                                     // ----------- NormalizeVector --------- //
  // ---------------------------------- //
  // Compute the angle between two vectors using the dot product and norms.
  // The angle is computed using the formula:
  // theta = acos((v1 . v2) / (||v1|| * ||v2||))
  // where v1 and v2 are the vectors, ||v1|| and ||v2|| are the norms of the vectors,
  // and theta is the angle between the vectors.
  // The angle is returned in radians.
  // ---------------------------------- //
  T VectorAngle(const Vec &v1, const Vec &v2) const
  {                                     // ----------- VectorAngle ------------ //
    if (v1.size()!=v2.size())           // Are they the same length?
    {                                   // No, can't compute angle.
      std::cerr<<"[ERROR] Vectors must have the same length for angle computation!\n";
      return T(0);                      // Return 0 if they are not the same length.
    }                                   // Done checking if they are the same length.
    T dot=VectorDotProduct(v1,v2);      // Compute the dot product of the vectors.
    T norm1=VectorNorm(v1);             // Compute the norm of the first vector.
    T norm2=VectorNorm(v2);             // Compute the norm of the second vector.
    if (norm1==T(0)||norm2==T(0))       // Are any of the norms zero?
    {                                   // Yes, can't compute angle.
      std::cerr<<"[ERROR] Cannot compute angle with zero norm vector!\n";
      return T(0);                      // Return 0 if any of the norms is zero.
    }                                   // Done checking if any of the norms is zero.
    return std::acos(dot/(norm1*norm2)); // Return the angle in radians.
  }                                     // ----------- VectorAngle ------------ //
  // ---------------------------------- //
  // Compute the cross product of two vectors. The cross product is only defined
  // for 3-dimensional vectors. It is a vector that is orthogonal to both input vectors.
  // The cross product is used in physics to compute the torque and angular momentum.
  // It is used in spectral analysis to compute the orthogonal basis of a set of vectors.
  // The cross product is defined as:
  // v1 x v2 = (v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0])
  // where v1 and v2 are the vectors.
  // The function returns a new vector that is the cross product of the input vectors.
  // ---------------------------------- //
  Vec VectorCrossProduct(const Vec &v1, const Vec &v2) const
  {                                     // ----------- VectorCrossProduct ------ //
    if (v1.size()!=3||v2.size()!=3)     // Are they both 3-dimensional?
    {                                   // No, can't compute cross product.
      std::cerr<<"[ERROR] Vectors must be 3-dimensional for cross product!\n";
      return Vec{};                     // Return an empty vector if they are not 3-dimensional.
    }                                   // Done checking if they are both 3-dimensional.
    Vec cross(3);                       // Create a new vector of size 3.
    cross[0]=v1[1]*v2[2]-v1[2]*v2[1];   // Compute the first element of the cross product.
    cross[1]=v1[2]*v2[0]-v1[0]*v2[2];   // Compute the second element of the cross product.
    cross[2]=v1[0]*v2[1]-v1[1]*v2[0];   // Compute the third element of the cross product.
    return cross;                       // Return the cross product vector.
  }                                     // ----------- VectorCrossProduct ------ //
  // ---------------------------------- //
  // Add two vectors element-wise.
  // The function returns a new vector that is the sum of the two input vectors.
  // The input vectors must be of the same size.
  // ----------------------------------- //
  Vec AddVectors(const Vec &v1, const Vec &v2) const
  {                                     // ----------- AddVectors ------------ //
    if (v1.size()!=v2.size())           // Are they the same length?
    {                                   // No, can't add vectors.
      std::cerr<<"[ERROR] Vectors must have the same length for addition!\n";
      return Vec{};                     // Return an empty vector if they are not the same length.
    }                                   // Done checking if they are the same length.
    Vec sum(v1.size());                 // Create a new vector of the same size as v1 and v2.
    for (size_t i=0;i<v1.size();i++)    // For each element in the vectors...
    {                                   // Loop over the elements of the vectors.
      sum[i]=v1[i]+v2[i];               // Add the elements of the vectors.
    }                                   // Done adding the vectors.
    return sum;                         // Return the sum of the vectors.
  }                                     // ----------- AddVectors ------------ //
  // ---------------------------------- //
  // Subtract two vectors element-wise.
  // The function returns a new vector that is the difference of the two input vectors.
  // The input vectors must be of the same size.
  // ----------------------------------- //
  Vec SubtractVectors(const Vec &v1, const Vec &v2) const
  {                                     // ----------- SubtractVectors ------- //
    if (v1.size()!=v2.size())           // Are they the same length?
    {                                   // No, can't subtract vectors.
      std::cerr<<"[ERROR] Vectors must have the same length for subtraction!\n";
      return Vec{};                     // Return an empty vector if they are not the same length.
    }                                   // Done checking if they are the same length.
    Vec diff(v1.size());                // Create a new vector of the same size as v1 and v2.
    for (size_t i=0;i<v1.size();i++)    // For each element in the vectors...
    {                                   // Loop over the elements of the vectors.
      diff[i]=v1[i]-v2[i];              // Subtract the elements of the vectors.
    }                                   // Done subtracting the vectors.
    return diff;                        // Return the difference of the vectors.
  }                                     // ----------- SubtractVectors ------- //
  // ---------------------------------- //
  // Scale a vector by a scalar value.
  // ---------------------------------- //S
  Vec ScaleVector(const Vec &v, T scalar) const
  {                                     // ----------- ScaleVector ------------ //
    if (v.empty())                      // Is the vector emtpy?
    {                                   // Yes, return an empty vector.
      std::cerr<<"[ERROR] Vector is emtpy!\n";
      return Vec{};                     // Return an empty vector if the vector is empty.
    }                                   // Done checking if the vector is empty.
    Vec scaled(v.size());               // Create a new vector of the same size as v.
    for (size_t i=0;i<v.size();i++)     // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      scaled[i]=v[i]*scalar;            // Scale each element by the scalar value.
    }                                   // Done scaling the vector.
    return scaled;                      // Return the scaled vector.
  }                                     // ----------- ScaleVector ------------ //

private:
  Vec vect;                             // The vector of type T.
  
};
#endif // VECTORS_H
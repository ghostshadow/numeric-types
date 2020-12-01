/// \file

#ifndef NUMERICTYPES_MAT_HPP
#define NUMERICTYPES_MAT_HPP

#include <type_traits>
#include <stdexcept>
#include <cmath>
#include "Vec.hpp"

namespace NumericTypes {

/* ********** Matrix ********* */

template<class T, size_t R, size_t C>
class Mat;

// some sfine magic ...
/**
 * \brief Determinant of square matrix
 * \note This function is implemented for:
 * * 0x0-Matrix (defined as 1.)
 * * 1x1-Matrix (determinat of a scalar is itself)
 * * 2x2-Matrix
 * * 3x3-Matrix
 * * NxN-Matrix (for all matricies larger then 3,
 *   implemented using laplace expansion)
 * \note The correct implementation is chosen at compile time via SFINAE std::enable_if
 * \tparam T Matrix elements type
 * \tparam N Matrix dimensions
 * \param mat The matrix
 * \return Determinant of matrix
 */
template<class T, size_t N, typename std::enable_if<(N>3), bool>::type=true>
inline T determinant(const Mat<T, N, N> &mat) {
	T det{0.};
	for(size_t i=0; i<N; ++i) {
		det+=(1-2*(i%2))*mat(0, i)*mat.subMatrix(0, i).determinant();
	}
	return det;
}

// documented in NxN case
template<class T, size_t N, typename std::enable_if<(N==0), bool>::type=true>
inline T determinant(const Mat<T, N, N> &mat) {
	return 1.;
}

// documented in NxN case
template<class T, size_t N, typename std::enable_if<(N==1), bool>::type=true>
inline T determinant(const Mat<T, N, N> &mat) {
	return mat(0, 0);
}

// documented in NxN case
template<class T, size_t N, typename std::enable_if<(N==2), bool>::type=true>
inline T determinant(const Mat<T, N, N> &mat) {
	return mat(0, 0)*mat(1, 1)-mat(0, 1)*mat(1, 0);
}

// documented in NxN case
template<class T, size_t N, typename std::enable_if<(N==3), bool>::type=true>
inline T determinant(const Mat<T, N, N> &mat) {
	return mat(0, 0)*mat(1, 1)*mat(2, 2)
		   +mat(0, 1)*mat(1, 2)*mat(2, 0)
		   +mat(0, 2)*mat(1, 0)*mat(2, 1)
		   -mat(0, 2)*mat(1, 1)*mat(2, 0)
		   -mat(0, 0)*mat(1, 2)*mat(2, 1)
		   -mat(0, 1)*mat(1, 0)*mat(2, 2);
}

/**
 * \brief Inverse of square matrix
 * \note This function is implemented for:
 * * 1x1-Matrix (inverse of a scalar is 1 divided by the scalar)
 * * NxN-Matrix (for all matricies larger then 1,
 *   implemented via the ajungate)
 * \note The correct implementation is chosen at compile time via SFINAE std::enable_if
 * \tparam T Matrix elements type
 * \tparam N Matrix dimensions
 * \param mat The matrix
 * \return A new Mat being the inverse of the matrix
 */
template<class T, size_t N, typename std::enable_if<N!=1, bool>::type=true>
inline Mat<T, N, N> inverse(const Mat<T, N, N> &mat) {
	static_assert(N>=1, "inverse does not make sense for zero-sized matrix");
	Mat<T, N, N> coef{};
	for(size_t r=0; r<N; ++r) {
		for(size_t c=0; c<N; ++c) {
			coef(r, c)=(1-2*(r%2))*(1-2*(c%2))*mat.subMatrix(r, c).determinant();
		}
	}
	return coef.transposed()/mat.determinant();
}

// documented in NxN case
template<class T, size_t N, typename std::enable_if<(N==1), bool>::type=true>
inline Mat<T, N, N> inverse(const Mat<T, N, N> &mat) {
	return Mat<T, 1, 1>{1./mat(0, 0)};
};

/**
 * \brief General RxC-dimensional Matrix
 * \tparam T Elements type
 * \tparam R Row count
 * \tparam C Column count
 * \pre The elements type has to support common arithmetic operations
 */
template<class T=double, size_t R=3, size_t C=3>
class Mat {
	static_assert(std::is_arithmetic<T>::value,
			"matrix elements must support arithmetic operations");
public:
	// Constructors
	/**	\brief Default constructor (all elements 0) */
	Mat()=default;

	/**
	 * \brief Construct Mat from array
	 * \details Fill the matrix row by row with the elements of the array
	 * \tparam T2 Array elements type
	 * \param _v The Array
	 * \pre 1. The array elements type has to be convertible to the matrix elements type
	 * \pre 2. The array has to have as many elements as the matrix does entries (R*C)
	 */
	template<class T2=T>
	explicit Mat(const T2 (&_v)[R*C]) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=_v[c+r*C];
			}
		}
	}

	/**
	 * \brief Construct Mat from initializer list
	 * \details Fill the matrix row by row with the elements of the initializer
	 * list. All none filled elements are set to 0, any excess elements are
	 * discarded
	 * \tparam T2 Initializer list elements type
	 * \param il The initializer list
	 * \pre The initializer list element type has to be convertible to the matrix
	 * elements type
	 */
	template<class T2=T>
	Mat(std::initializer_list<T2> il) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(int i=0; i<(il.size()>=(R*C)?(R*C):il.size()); ++i) {
			v[i]=*(il.begin()+i);
		}
	}

	/**
	 * \brief Construct Mat from array of Vec
	 * \details Fills the columns of the matrix with the vectors in the array
	 * \tparam T2 Vec elements type
	 * \param cols Column vectors array
	 * \pre 1. The elements type of the Vecs have to be convertible to the matrix
	 * elements type
	 * \pre 2. The Vecs have to be as large as the matrix has rows and the
	 * array has to habe as many elements as the matrix columns
	 */
	template<class T2=T>
	explicit Mat(const Vec<T2, R> (&cols)[C]) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=cols[c][r];
			}
		}
	}

	/**
	 * \brief Construct Mat from other Mat with different elements type
	 * \tparam T2 Elements type of the other Mat
	 * \param other The other Mat
	 * \pre 1. The elements type of the other Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The other Mat has to have the same dimensions as the matrix
	 * (identical row and column count)
	 */
	template<class T2=T>
	explicit Mat(const Mat<T2, R, C> &other) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=other(r, c);
			}
		}
	}

	/**
	 * \brief Construct Mat copy
	 * \param other The other Mat instance
	 */
	Mat(const Mat<T, R, C> &other) noexcept {
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=other(r, c);
			}
		}
	}

	// Predefinded Matricies
	/**
	 * \brief Create a identity matrix
	 * \details Creates a matrix with all elements 0 except the main diagonal,
	 * which is all 1s
	 * \return Identity matrix
	 */
	static Mat<T, R, C> Identity() noexcept {
		Mat<T, R, C> diag{};
		for(int i=0; i<(R>=C?C:R); ++i) {
			diag(i, i)=1.;
		}
		return diag;
	}

	/**
	 * \brief Create a diagonal matrix from Vec
	 * \details Creates a matrix with all elements 0 and the main diagonal
	 * (where both indices are equal) is filled with the elements of the
	 * Vec
	 * \tparam T2 Elements type of the Vec
	 * \tparam N Number of elements in the Vec (its dimension)
	 * \param vec The Vec
	 * \return The diagonal matrix
	 * \pre 1. The elements type of the Vec has to be convertible to the matrix
	 * elements type
	 * \pre 2. The Vec has to have at least as many elements as the diagonal
	 * has entries
	 */
	template<class T2=T, size_t N=3>
	static Mat<T, R, C> Diagonal(const Vec <T2, N> &vec) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		static_assert(N>=(R>=C?C:R),
				"vector for filling the diagonal has to be as long as the diagonal");
		Mat<T, R, C> diag{};
		for(int i=0; i<(R>=C?C:R); ++i) {
			diag(i, i)=vec[i];
		}
		return diag;
	}

	/**
	 * \brief Create a Mat representing a Vec as column vector
	 * \tparam T2 Elements type of the Vec
	 * \tparam N Number of elements of the Vec (its dimension)
	 * \param vec The Vec
	 * \return A matrix with N rows and only one column
	 * \pre The elements type of the Vec has to be convertible to the matrix
	 * elements type
	 */
	template<class T2=T, size_t N=3>
	static Mat<T, N, 1> ColumnVector(const Vec <T2, N> &vec) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, N, 1> col{};
		for(int i=0; i<N; ++i) {
			col(i, 0)=vec[i];
		}
		return col;
	}

	/**
	 * \brief Create a Mat representing a Vec as row vector
	 * \tparam T2 Elements type of the Vec
	 * \tparam N Number of elements of the Vec (its dimension)
	 * \param vec The Vec
	 * \return A matrix with only one row and N columns
	 * \pre The elements type of the Vec has to be convertible to the matrix
	 * elements type
	 */
	template<class T2=T, size_t N=3>
	static Mat<T, 1, N> RowVector(const Vec <T2, N> &vec) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, 1, N> col{};
		for(int i=0; i<N; ++i) {
			col(0, i)=vec[i];
		}
		return col;
	}

	// Index Accessors
	/**
	 * \brief Index access for modifyable lvalue
	 * \param row Row index (0-(nr of rows-1))
	 * \param column Column index (0-(nr of columns-1))
	 * \return Modifyable lvalue reference of the element at indices
	 * \throws std::out_of_range if any of the two indices are outside the
	 * bounds
	 */
	T &operator()(const size_t row, const size_t column) {
		if(row>=R || column>=C) {
			throw std::out_of_range("index out of bound");
		}
		return v[column+row*C];
	}

	/**
	 * \brief Index access for copy of element
	 * \param row Row index (0-(nr of rows-1))
	 * \param column Column index (0-(nr of columns-1))
	 * \return Copy of the element at indices
	 * \throws std::out_of_range if any of the two indices are outside the
	 * bounds
	 */
	T operator()(const size_t row, const size_t column) const {
		if(row>=R || column>=C) {
			throw std::out_of_range("index out of bound");
		}
		return v[column+row*C];
	}

	// Modifying Operators
	/**
	 * \brief Asignment from array
	 * \details Fill the matrix row by row with the elements of the array
	 * \tparam T2 Elements type of the array
	 * \param _v The array
	 * \return lvalue reference of self
	 * \pre 1. The array elements type has to be convertible to the matrix elements type
	 * \pre 2. The array has to have as many elements as the matrix does entries (R*C)
	 */
	template<class T2=T>
	Mat<T, R, C> &operator=(const T2 (&_v)[R*C]) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=_v[c+r*C];
			}
		}
		return *this;
	}

	/**
	 * \brief Asign from other Mat with different elements type
	 * \tparam T2 Elements type of the other Mat
	 * \param other The other Mat
	 * \return lvalue reference of self
	 * \pre 1. The elements type of the other Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The other Mat has to have the same dimensions as the matrix
	 * (identical row and column count)
	 */
	template<class T2=T>
	Mat<T, R, C> &operator=(const Mat<T2, R, C> &other) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=other(r, c);
			}
		}
		return *this;
	}

	/**
	 * \brief Asign copy
	 * \param other The other Mat instance
	 * \return lvalue reference of self
	 */
	Mat<T, R, C> &operator=(const Mat<T, R, C> &other) noexcept {
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]=other(r, c);
			}
		}
		return *this;
	}

	/**
	 * \brief Short asign sum
	 * \details Asigns each element the sum of itself and the corresponding
	 * element of the rhs Mat
	 * \tparam T2 Elements type of the rhs Mat
	 * \param rhs The rhs Mat
	 * \return lvalue reference of self
	 * \pre 1. The elements type of the rhs Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The rhs Mat has to have the same dimensions as the matrix
	 * (identical row and column count)
	 */
	template<class T2=T>
	Mat<T, R, C> &operator+=(const Mat<T2, R, C> &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]+=rhs(r, c);
			}
		}
		return *this;
	}

	/**
	 * \brief Short asign difference
	 * \details Asigns each element the difference of itself and the corresponding
	 * element of the rhs Mat
	 * \tparam T2 Elements type of the rhs Mat
	 * \param rhs The rhs Mat
	 * \return lvalue reference of self
	 * \pre 1. The elements type of the rhs Mat has to be convertible to the
	 * matrix elements typ
	 * \pre 2. The rhs Mat has to have the same dimensions as the matrix
	 * (identical row and column count)e
	 */
	template<class T2=T>
	Mat<T, R, C> &operator-=(const Mat<T2, R, C> &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]-=rhs(r, c);
			}
		}
		return *this;
	}

	/**
	 * \brief Short asign multiplicative scale
	 * \details Asigns each element the product of itself and the scaling value
	 * \tparam T2 The scaling value type
	 * \param rhs The scaling value
	 * \return lvalue reference of self
	 * \pre The scaling value type has to be convertible to the matrix
	 * elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
												 bool>::type=true>
	Mat<T, R, C> &operator*=(const T2 &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]*=rhs;
			}
		}
		return *this;
	}

	/**
	 * \brief Short asign dividing scale
	 * \details Asigns each element the product of itself and inverse of the scaling value
	 * \tparam T2 The scaling value type
	 * \param rhs The scaling value
	 * \return lvalue reference of self
	 * \pre The scaling value type has to be convertible to the matrix
	 * elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
												 bool>::type=true>
	Mat<T, R, C> &operator/=(const T2 &rhs) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				v[c+r*C]/=rhs;
			}
		}
		return *this;
	}

	// Non-Modifying Operators
	/**
	 * \brief Equality test
	 * \tparam T2 Elements type of the rhs Mat
	 * \param rhs The rhs Mat
	 * \return true if all elements of the matrix are equal to the corresponding
	 * elements of the rhs Mat
	 * \pre 1. The elements type of the rhs Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The rhs Mat has to have the same dimensions as the matrix
	 * (identical row and column count)
	 */
	template<class T2=T>
	bool operator==(const Mat<T2, R, C> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(int r=0; r<R; ++r) {
			for(int c=0; c<C; ++c) {
				if(v[c+r*C]!=rhs(r, c)) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * \brief Unary Negative
	 * \return A new Mat instance with each element being the negative of the corresponding
	 * element of the matrix
	 */
	Mat<T, R, C> operator-() const noexcept {
		Mat<T, R, C> neg{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				neg(r, c)=-v[c+r*C];
			}
		}
		return neg;
	}

	/**
	 * \brief Summation
	 * \tparam T2 Elements type of the rhs Mat
	 * \param rhs The rhs Mat
	 * \return A new Mat instance with each element being the sum of the
	 * corresponding elements of the matrix and the rhs Mat
	 * \pre 1. The elements type of the rhs Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The rhs Mat has to have the same dimensions as the matrix
	 * (identical row and column count)
	 */
	template<class T2=T>
	Mat<T, R, C> operator+(const Mat<T2, R, C> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, R, C> sum{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				sum(r, c)=v[c+r*C]+rhs(r, c);
			}
		}
		return sum;
	}

	/**
	 * \brief Difference
	 * \tparam T2 Elements type of the rhs Mat
	 * \param rhs The rhs Mat
	 * \return A new Mat instance with each element being the difference of the
	 * corresponding elements of the matrix and the rhs Mat
	 * \pre 1. The elements type of the rhs Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The rhs Mat has to have the same dimensions as the matrix
	 * (identical row and column count)
	 */
	template<class T2=T>
	Mat<T, R, C> operator-(const Mat<T2, R, C> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, R, C> diff{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				diff(r, c)=v[c+r*C]-rhs(r, c);
			}
		}
		return diff;
	}

	/**
	 * \brief Multiplicative scaled
	 * \tparam T2 The scaling value type
	 * \param rhs The scaling value
	 * \return A new Mat instance with each element being the product of the
	 * corresponding element of the matrix and the scaling value
	 * \pre The scaling value type has to be a scalar, which is convertible to the
	 * matrix elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
												 bool>::type=true>
	Mat<T, R, C> operator*(const T2 &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, R, C> scaled{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				scaled(r, c)=v[c+r*C]*rhs;
			}
		}
		return scaled;
	}

	/**
	 * \brief Product
	 * \details The matrix product is defined such that each element of the
	 * result matrix is the sum of the product of the corresponding elements
	 * of the row of the matrix and the column of the rhs Mat.
	 * \tparam T2 Elements type of the rhs Mat
	 * \tparam C3 Column count of the rhs Mat
	 * \param rhs The rhs Mat
	 * \return A new Mat instance with the row count being the row count of the
	 * matrix and the column count being the column count of the rhs Mat.
	 * The elements are set according to the matrix multiplication scheme.
	 * \pre 1. The elements type of the rhs Mat has to be convertible to the
	 * matrix elements type
	 * \pre 2. The rhs Mat has to have as many rows as the matrix has columns
	 */
	template<class T2=T, size_t C2=3>
	Mat<T, R, C2> operator*(const Mat<T2, C, C2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, R, C2> prod{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C2; ++c) {
				prod(r, c)=0.;
				for(size_t i=0; i<C; ++i) {
					prod(r, c)+=v[i+r*C]*rhs(i, c);
				}
			}
		}
		return prod;
	}

	/**
	 * \brief Dividing scaled
	 * \tparam T2 The scaling value type
	 * \param rhs The scaling value
	 * \return A new Mat instance with each element being the product of the
	 * corresponding element of the matrix and inverse of the scaling value
	 * \pre The scaling value type has to be a scalar, which is convertible to the
	 * matrix elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
												 bool>::type=true>
	Mat<T, R, C> operator/(const T2 &rhs) const {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Mat<T, R, C> scaled{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				scaled(r, c)=v[c+r*C]/rhs;
			}
		}
		return scaled;
	}

	// Named Operators
	/**
	 * \brief The Transposed
	 * \return A new Mat instance with dimensions and conent flipped along the
	 * diagonal
	 */
	Mat<T, C, R> transposed() const noexcept {
		Mat<T, C, R> transp{};
		for(size_t r=0; r<R; ++r) {
			for(size_t c=0; c<C; ++c) {
				transp(r, c)=v[c+r*C];
			}
		}
		return transp;
	}

	/**
	 * \brief The sub matrix
	 * \param rrow Index of the row to be removed
	 * \param rcolumn Index of the column to be removed
	 * \return A new Mat instance with dimensions 1 less than the matrix, filled
	 * such that the values of the specified row and colum have been removed and replaced
	 * with the elements of the next higher index
	 * \throws std::out_of_range if any of the two indices are outside the
	 * bounds
	 * \pre The matrix has to have more than one row and column
	 */
	Mat<T, R-1, C-1> subMatrix(const size_t rrow, const size_t rcolumn) const {
		static_assert(R>1 && C>1, "Can not create submatrix of Scalar value");
		if(rrow>=R || rcolumn>=C) {
			throw std::out_of_range("index out of bounds");
		}
		Mat<T, R-1, C-1> sub{};
		size_t ro=0, co=0;
		for(size_t r=0; r<R-1; ++r) {
			co=0;
			for(size_t c=0; c<C-1; ++c) {
				sub(r, c)=v[co+ro*C];
				co=(co+1)==rcolumn?co+2:co+1;
			}
			ro=(ro+1)==rrow?ro+2:ro+1;
		}
		return sub;
	}

	/**
	 * \brief The sub matrix (remove only a row)
	 * \param rrow Index of the row to be removed
	 * \return A new Mat instance with one less rows and the same number of columns
	 * than the matrix, filled such that the values of the specified row have been
	 * removed and replaced with the elements of the next higher index
	 * \throws std::out_of_range if the index is outside the bounds
	 * \pre The matrix has to have more than one row
	 */
	Mat<T, R-1, C> subMatrixRemoveRow(const size_t rrow) const {
		static_assert(R>1, "Can not remove Row, there is only one");
		if(rrow>=R) {
			throw std::out_of_range("index out of bounds");
		}
		Mat<T, R-1, C> sub{};
		size_t ro=0;
		for(size_t r=0; r<R-1; ++r) {
			for(size_t c=0; c<C; ++c) {
				sub(r, c)=v[c+ro*C];
			}
			ro=(ro+1)==rrow?ro+2:ro+1;
		}
		return sub;
	}

	/**
	 * \brief The sub matrix
	 * \param rcolumn Index of the column to be removed
	 * \return A new Mat instance with one less column than the matrix, filled
	 * such that the values of the specified colum have been removed and replaced
	 * with the elements of the next higher index
	 * \throws std::out_of_range if the index is outside the bounds
	 * \pre The matrix has to have more than one column
	 */
	Mat<T, R, C-1> subMatrixRemoveColumn(const size_t rcolumn) const {
		static_assert(C>1, "Can not remove Column, there is only one");
		if(rcolumn>=C) {
			throw std::out_of_range("index out of bounds");
		}
		Mat<T, R, C-1> sub{};
		size_t co=0;
		for(size_t r=0; r<R; ++r) {
			co=0;
			for(size_t c=0; c<C-1; ++c) {
				sub(r, c)=v[co+r*C];
				co=(co+1)==rcolumn?co+2:co+1;
			}
		}
		return sub;
	}

	/**
	 * \brief Get a copy of the specified column as Vec
	 * \param column Index of the column to get
	 * \return A new Vec instance containing the elements of the requested
	 * column
	 * \throws std::out_of_range if the index is out of bound
	 */
	Vec <T, R> getColumn(const size_t column) const {
		if(column>=C) {
			throw std::out_of_range("index out of bounds");
		}
		Vec<T, R> col{};
		for(int i=0; i<R; ++i) {
			col[i]=v[column+i*C];
		}
		return col;
	}

	/**
	 * \brief Get a copy of the specified row as Vec
	 * \param row Index of the row to get
	 * \return A new Vec instance containing the elements of the requested
	 * row
	 * \throws std::out_of_range if the index is out of bound
	 */
	Vec <T, C> getRow(const size_t row) const {
		if(row>=R) {
			throw std::out_of_range("index out of bounds");
		}
		Vec<T, C> rowvec{};
		for(int i=0; i<C; ++i) {
			rowvec[i]=v[i+row*C];
		}
		return rowvec;
	}

	// Square Matrix Named Operators
	/**
	 * \brief Trace of the matrix
	 * \return Sum of the elements in the main diagonal
	 * \pre The matrix has to have as many columns as rows
	 */
	T trace() const noexcept {
		static_assert(R==C, "trace is only defined for square matricies");
		T trace{0.};
		for(size_t i=0; i<C; ++i) {
			trace+=v[i+i*C];
		}
		return trace;
	}

	/**
	 * \brief Determinant of the matrix
	 * \return Determinant according to the determinant definition of a matrix
	 * \pre The matrix has to have as many columns as rows
	 */
	T determinant() const noexcept {
		static_assert(R==C, "determinant is only defined for square matricies");
		return NumericTypes::determinant(*this);
	}

	/**
	 * \brief Inverse of the matrix
	 * \details The inverse is a matrix, which if multiplied with the matrix
	 * results in the identity matrix
	 * \return The inverse of the matrix according to the inverse definition of
	 * a matrix
	 * \pre The matrix has to have as many columns as rows
	 */
	Mat<T, R, C> inverse() const {
		static_assert(R==C, "inverse is only defined for square matricies");
		return NumericTypes::inverse(*this);
	}

	// Type Info
	/** \brief Row count */
	const size_t rows=R;
	/** \brief Column count */
	const size_t columns=C;
	/** \brief Elements type */
	using element_type=T;

private:
	/** \brief instance elements storage */
	T v[R*C]{0.};
};

/**
 * \brief Multiplicative scaled
 * \tparam T Elements type
 * \tparam R Row count
 * \tparam C Column count
 * \tparam T2 The scaling value type
 * \param lhs The scaling value
 * \param rhs The matrix
 * \return A new Mat instance with each element being the product of the
 * corresponding element of the matrix and the scaling value
 * \pre The scaling value type has to be a scalar, which is convertible to the
 * matrix elements type
 */
template<class T=double,
		 class T2=T, size_t R=3, size_t C=3, typename std::enable_if<std::is_scalar<T2>::value,
																	 bool>::type=true>
inline Mat<T, R, C> operator*(const T2 &lhs, const Mat<T, R, C> &rhs) noexcept {
	static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	return rhs*lhs;
}

/**
 * \brief Automatic Product Mat-Vec
 * \tparam T Vec elements type
 * \tparam T2 Mat elements type
 * \tparam R Mat row count
 * \tparam C Mat column count
 * \param lhs The lhs Mat
 * \param rhs The rhs Vec (interpreted as column vector)
 * \return A new Vec instance with as many elements as the Mat has rows and each
 * element being the sum of the product of the corresponding lhs Mat rows and
 * the Vec element
 * \pre The Mat elements type has to be convertible to the Vec elements type
 * \pre The Vec has to have as many elements as the Mat has columns
 */
template<class T=double, class T2=T, size_t R=3, size_t C=3>
inline Vec <T, R> operator*(const Mat<T2, R, C> &lhs, const Vec <T, C> &rhs) noexcept {
	static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	return (lhs*Mat<T2, R, C>::ColumnVector(rhs)).getColumn(0);
};

/**
 * \brief Automatic Product Vec-Mat
 * \tparam T Vec elements type
 * \tparam T2 Mat elements type
 * \tparam R Mat row count
 * \tparam C Mat column count
 * \param lhs The lhs Vec
 * \param rhs The rhs Mat
 * \return A new Vec instance with as many elements as the Mat has columns and each
 * element being the sum of the product of the corresponding lhs Mat columns and
 * the Vec element
 * \pre The Mat elements type has to be convertible to the Vec elements type
 * \pre The Vec has to have as many elements as the Mat has rows
 */
template<class T=double, class T2=T, size_t R=3, size_t C=3>
inline Vec <T, C> operator*(const Vec <T, R> &lhs, const Mat<T2, R, C> &rhs) noexcept {
	static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	return (Mat<T2, R, C>::RowVector(lhs)*rhs).getRow(0);
};

/** \brief Shorthand for 3x3-Matrix of type double */
using Mat33D=Mat<double, 3, 3>;
/** \brief Shorthand for 4x4-Matrix of type double */
using Mat44D=Mat<double, 4, 4>;

}

#endif //NUMERICTYPES_MAT_HPP

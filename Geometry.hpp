/// \file

#ifndef NUMERICTYPES_GEOMETRY_HPP
#define NUMERICTYPES_GEOMETRY_HPP

#include "Vec.hpp"
#include "Mat.hpp"
#include "Quat.hpp"

namespace NumericTypes {

/**
 * \name Geometry Helper Functions
 * \brief Functions to help with 3D-geometric operations
 * @{
 */

/**
 * \brief Rotation axis for 3D-rotations
 */
enum class RotationAxis {
	X /**< \brief Rotation around X-axis */,
	Y /**< \brief Rotation around Y-axis */,
	Z /**< \brief Rotation around Z-axis */
};

/**
 * \brief Create a rotation matrix for a 3D-rotation around a single axis
 * \tparam T Mat elements type
 * \tparam T2 Rotation angle value type
 * \param ax Axis around which to rotate (see RotationAxis)
 * \param ang_rad Angle to rotate (sign according to right hand rule) [(rad)]
 * \return 3D-rotation matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \pre The angle value type has to be convertible to the Mat elements type
 */
template<class T=double, class T2=T,
	typename std::enable_if<std::is_scalar<T2>::value, bool>::type=true>
inline Mat<T, 3, 3> RotationMatrix(const RotationAxis ax, const T2 ang_rad) {
	static_assert(std::is_convertible<T2, T>::value,
			"angle type has to be convertible to matrix element type");
	switch(ax) {
	case RotationAxis::X:
		return Mat<T, 3, 3>{
			1., 0., 0.,
			0., std::cos(ang_rad), std::sin(ang_rad),
			0., -sin(ang_rad), cos(ang_rad)
		};
	case RotationAxis::Y:
		return Mat<T, 3, 3>{
			std::cos(ang_rad), 0., -std::sin(ang_rad),
			0., 1., 0.,
			std::sin(ang_rad), 0., std::cos(ang_rad)
		};
	case RotationAxis::Z:
		return Mat<T, 3, 3>{
			std::cos(ang_rad), std::sin(ang_rad), 0.,
			-std::sin(ang_rad), std::cos(ang_rad), 0.,
			0., 0., 1.
		};
	}
}

/**
 * \brief Rotate a vector to the new attitude represented by quaternion
 * \tparam T Vec elements type
 * \tparam T2 Quat elements type
 * \param v The Vec to be rotated
 * \param q The Quat representing the rotation (should be a versor)
 * \return Rotated vector
 * \pre The Quat elements type has to be convertible to the Vec elements type
 */
template<class T=double, class T2=T>
inline Vec<T, 3> rotateVecByQuat(const Vec<T, 3> &v, const Quat <T2> &q) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	return (q.normalized()*Quat<T>(v)*q.normalized().reciprocal()).vector();
}

/**
 * \brief Create a quaternion representing a rotation around an arbitrary axis
 * \tparam T Quat elements type
 * \tparam T2 Axis Vec elements type
 * \tparam T3 Angle value type
 * \param a Vec representing the axis (also determines the positive rotation direction
 * via the right hand rule)
 * \param theta Angle to rotate [(rad)]
 * \return The Quat representing the rotation
 * \pre The axis Vec elements type and the angel value type has to be convertible
 * to the Quat elements type
 */
template<class T=double, class T2=T, class T3=T,
	typename std::enable_if<std::is_scalar<T3>::value, bool>::type=true>
inline Quat<T> quatFromAxisAndRoatation(const Vec<T2, 3> &a, const T3 theta) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	static_assert(std::is_convertible<T3, T>::value, "element type is not convertible");
	if(a.norm2() == 0)
		throw std::domain_error("Invalid rotation axis");
	return Quat<T>(cos(theta/2), a.normalized()*(sin(theta/2))).normalized();
}

/**
 * \brief Convert attitude from quaternion to matrix representation
 * \tparam T Mat elements type
 * \tparam T2 Quat elements type
 * \param quat Quaternion representing the attitude (should be a versor)
 * \return The attitude matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \pre The Quat elements type has to be convertible to the Mat elements type
 */
template<class T=double, class T2=T>
inline Mat<T, 3, 3> attitudeMatrixFromQuat(const Quat <T2> &quat) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	Quat<T2> nq=quat.normalized();
	using Q=QuatElement;
	return Mat<T, 3, 3>{
			nq[Q::S]*nq[Q::S]+nq[Q::X]*nq[Q::X]-nq[Q::Y]*nq[Q::Y]-
				nq[Q::Z]*nq[Q::Z], /*00*/
			2.*(nq[Q::X]*nq[Q::Y]+nq[Q::S]*nq[Q::Z]), /*01*/
			2.*(nq[Q::X]*nq[Q::Z]-nq[Q::S]*nq[Q::Y]), /*02*/
			2.*(nq[Q::X]*nq[Q::Y]-nq[Q::S]*nq[Q::Z]), /*10*/
			nq[Q::S]*nq[Q::S]-nq[Q::X]*nq[Q::X]+nq[Q::Y]*nq[Q::Y]-
				nq[Q::Z]*nq[Q::Z], /*11*/
			2.*(nq[Q::Y]*nq[Q::Z]+nq[Q::S]*nq[Q::X]), /*12*/
			2.*(nq[Q::X]*nq[Q::Z]+nq[Q::S]*nq[Q::Y]), /*20*/
			2.*(nq[Q::Y]*nq[Q::Z]-nq[Q::S]*nq[Q::X]), /*21*/
			nq[Q::S]*nq[Q::S]-nq[Q::X]*nq[Q::X]-nq[Q::Y]*nq[Q::Y]+nq[Q::Z]*nq[Q::Z] /*22*/
	};
}

/**
 * \brief Convert attitude from matrix to quaternion represenatation
 * \tparam T Quat elements type
 * \tparam T2 Mat elements type
 * \param mat Attitude matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \return The attitude quaternion (is a versor)
 * \pre The Mat elements type has to be convertible to the Quat elements type
 */
template<class T, class T2=T>
inline Quat <T> quatFromAttitudeMatrix(const Mat<T2, 3, 3> &mat) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	T S{.5*sqrt(1.+(mat(0, 0)+mat(1, 1)+mat(2, 2)))};
	Vec<T, 3> V{
		(1./(4.*S))*(mat(1, 2)-mat(2, 1)), /*X*/
		(1./(4.*S))*(mat(2, 0)-mat(0, 2)), /*Y*/
		(1./(4.*S))*(mat(0, 1)-mat(1, 0)) /*Z*/
	};
	return Quat<T>(S, V).normalized();
}

/**
 * \brief Convert attitude from 321-Euler angles (z-y'-x") to matrix representation
 * \tparam T Mat elements type
 * \tparam T2 Vec elements type
 * \param e321 Vector containing the angles [(deg)] in turning order
 * (i.e. element 0 contains the angle to be turned around z, element 1 the angle
 * for y' and element 2 the angle for x")
 * \note It is up to the user to make sure the angles are not effected by
 * gimbal lock.
 * \return The attitude matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \pre The Vec elements type has to be convertible to the Mat elements type
 */
template<class T=double, class T2=T>
inline Mat<T, 3, 3> attitudeMatrixFromEuler321(const Vec<T2, 3> &e321) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	using V=VecElement;
	using RA=RotationAxis;
	// R(a1,a2,a3)=R(x",a3)*R(y',a2)*R(z,a1)
	// S=R(a1,a2,a3)*E
	return RotationMatrix(RA::X, M_PI/180.*e321[V::Z])*
		   RotationMatrix(RA::Y, M_PI/180.*e321[V::Y])*
		   RotationMatrix(RA::Z, M_PI/180.*e321[V::X]);
}

/**
 * \brief Convert attitude from matrix to 321-Euler angles (z-y'-x") representation
 * \tparam T Vec elements type
 * \tparam T2 Mat elements type
 * \param am The attitude matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \return A vector containing the angles [(deg)] in turning order
 * (i.e. element 0 is the angle around z, element 1 the angle around y' and
 * element 2 the angle around x")
 * \pre The Mat elements type has to be convertible to the Vec elements type
 */
template<class T=double, class T2=T>
inline Vec<T, 3> euler321FromAttitudeMatrix(const Mat<T2, 3, 3> &am) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	return (180/M_PI)*Vec<T, 3>{std::atan2(am(0, 1), am(0, 0)),
		-std::asin(am(0, 2)),
		std::atan2(am(1, 2), am(2, 2))};
}

/**
 * \brief Convert attitude from 313-Euler angles (z-x'-z") to matrix representation
 * \tparam T Mat elements type
 * \tparam T2 Vec elements type
 * \param e313 Vector containing the angles [(deg)] in turning order
 * (i.e. element 0 contains the angle to be turned around z, element 1 the angle
 * for x' and element 2 the angle for z")
 * \note It is up to the user to make sure the angles are not effected by
 * gimbal lock.
 * \return The attitude matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \pre The Vec elements type has to be convertible to the Mat elements type
 */
template<class T=double, class T2=T>
inline Mat<T, 3, 3> attitudeMatrixFromEuler313(const Vec<T2, 3> &e313) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	using V=VecElement;
	using RA=RotationAxis;
	// R(a1,a2,a3)=R(z",a3)*R(x',a2)*R(z,a1)
	// S=R(a1,a2,a3)*E
	return RotationMatrix(RA::Z, M_PI/180.*e313[V::Z])*
		   RotationMatrix(RA::X, M_PI/180.*e313[V::Y])*
		   RotationMatrix(RA::Z, M_PI/180.*e313[V::X]);
}

/**
 * \brief Convert attitude from matrix to 313-Euler angles (z-x'-z") representation
 * \tparam T Vec elements type
 * \tparam T2 Mat elements type
 * \param am The attitude matrix
 * (as the vector base of the new reference in terms of the previous reference as
 * line vectors)
 * \return A vector containing the angles [(deg)] in turning order
 * (i.e. element 0 is the angle around z, element 1 the angle around x' and
 * element 2 the angle around z")
 * \pre The Mat elements type has to be convertible to the Vec elements type
 */
template<class T=double, class T2=T>
inline Vec<T, 3> euler313FromAttitudeMatrix(const Mat<T2, 3, 3> &am) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	return (180/M_PI)*Vec<T, 3>{std::atan2(am(2, 0), -am(2, 1)),
		std::acos(am(2, 2)),
		std::atan2(am(0, 2), am(1, 2))};
}


/**
 * \brief Spherical linear interpolation of quaternions
 * \tparam T start Quat elements type
 * \tparam T2 end Quat elements type
 * \tparam T3 fraction type
 * \param start Interpolation range start quaternion
 * \param end Interpolation range end quaternion
 * \param t Interpolation fraction
 * (values between 0 and 1 lie between the two quaternions)
 * \return Interpolated quaternion
 */
template<class T=double, class T2=T, class T3=T,
	typename std::enable_if<std::is_scalar<T3>::value, bool>::type=true>
inline Quat<T> slerp(const Quat<T> &start, const Quat<T2> &end, const T3 t) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	static_assert(std::is_convertible<T3, T>::value, "element type is not convertible");
	T angle = std::acos(start.normalized().dot(end.normalized()));
	return Quat<T>((std::sin((1.-t)*angle)/std::sin(angle)) * start +
			(std::sin(t*angle)/std::sin(angle)) * end);
}

/**
 * \brief Linear interpolation of vectors
 * \tparam N vector dimensionality
 * (can only interpolate between vectors with same dimensions)
 * \tparam T start Vec elements type
 * \tparam T2 end Vec elements type
 * \tparam T3 fraction type
 * \param start Interpolation range start vector
 * \param end Interpolation range end vector
 * \param t Interpolation fraction
 * (values between 0 and 1 lie between the two vectors)
 * \return Interpolated vector
 */
template<size_t N=3, class T=double, class T2=T, class T3=T,
	typename std::enable_if<std::is_scalar<T3>::value, bool>::type=true>
inline Vec<T, N> lerp(const Vec<T, N> &start, const Vec<T2, N> &end, const T3 t) {
	static_assert(std::is_convertible<T2, T>::value, "element type is not convertible");
	static_assert(std::is_convertible<T3, T>::value, "element type is not convertible");
	return Vec<T, N>((1-t) * start + t * end);
}



/**
 * @}
 */

}


#endif //NUMERICTYPES_GEOMETRY_HPP

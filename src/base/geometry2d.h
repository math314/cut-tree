#pragma once
#include <complex>

namespace agl {
namespace geometry2d {
const double kEps = 1e-10;
using point_type = std::complex<double>;
using segment_type = std::pair<point_type, point_type>;

inline double dot(const point_type &p, const point_type &q) { return real(conj(p) * q); }
inline double det(const point_type &p, const point_type &q) { return imag(conj(p) * q); }

inline int ccw(const point_type &a, point_type b, point_type c) {
  b -= a; c -= a;
  if (det(b, c) > kEps) return 1;          // ccw
  if (det(b, c) < -kEps) return -1;        // cw
  if (dot(b, c) < -kEps) return 2;         // c-a-b
  if (norm(b) < norm(c) - kEps) return -2; // a-b-c
  return 0;                               // a-c-b, b==c
}

inline bool does_intersect(const segment_type &s1, const segment_type &s2) {
  return
      ccw(s1.first, s1.second, s2.first) * ccw(s1.first, s1.second, s2.second) <= 0 &&
      ccw(s2.first, s2.second, s1.first) * ccw(s2.first, s2.second, s1.second) <= 0;
}
}
}

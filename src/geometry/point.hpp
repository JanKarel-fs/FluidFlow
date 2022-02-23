#ifndef POINT_HPP
#define POINT_HPP

template <typename T>
class Point2 {
public:
  T x;
  T y;

  Point2() {};

  Point2(T X, T Y): x(X), y(Y) {};

  ~Point2() {};
};

template <typename T>
inline Point2<T> operator+(const Point2<T>& a, const Point2<T>& b) {
  return Point2<T>(a.x+b.x, a.y+b.y);
}

template <typename T, typename S>
inline Point2<T> operator*(const Point2<T>& a, const S& b) {
  return Point2<T>(a.x*b, a.y*b);
}

template <typename T, typename S>
inline Point2<T> operator*(const S& b, const Point2<T>& a) {
  return Point2<T>(a.x*b, a.y*b);
}

template <typename T, typename S>
inline Point2<T> operator/(const Point2<T>& a, const S& b) {
  return Point2<T>(a.x/b, a.y/b);
}

typedef Point2<double> Point2d;

#endif

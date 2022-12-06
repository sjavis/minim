#ifndef VEC_HPP
#define VEC_HPP

#include <math.h>
#include <cassert>
#include <numeric>
#include <algorithm>


// Sum
template<typename T, typename U>
auto operator+(T a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(b.size());
  std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a+x; });
  return c;
}
template<typename T, typename U>
auto operator+(const vector<T>& a, U b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x+b; });
  return c;
}
template<typename T, typename U>
auto operator+(const vector<T>& a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<T>());
  return c;
}

template<typename T, typename U>
auto operator+=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x+b; });
  return a;
}
template<typename T, typename U>
auto operator+=(vector<T>& a, const vector<U>& b) {
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<T>());
  return a;
}


// Diff
template<typename T>
auto operator-(const vector<T>& a) {
  vector<T> c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [](auto x){ return -x; });
  return c;
}

template<typename T, typename U>
auto operator-(T a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(b.size());
  std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a-x; });
  return c;
}
template<typename T, typename U>
auto operator-(const vector<T>& a, U b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x-b; });
  return c;
}
template<typename T, typename U>
auto operator-(const vector<T>& a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<T>());
  return c;
}

template<typename T, typename U>
auto operator-=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x-b; });
  return a;
}
template<typename T, typename U>
auto operator-=(vector<T>& a, const vector<U>& b) {
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<T>());
  return a;
}


// Multiply
template<typename T, typename U>
auto operator*(T a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(b.size());
  std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a*x; });
  return c;
}
template<typename T, typename U>
auto operator*(const vector<T>& a, U b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x*b; });
  return c;
}
template<typename T, typename U>
auto operator*(const vector<T>& a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<T>());
  return c;
}

template<typename T, typename U>
auto operator*=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x*b; });
  return a;
}
template<typename T, typename U>
auto operator*=(vector<T>& a, const vector<U>& b) {
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<T>());
  return a;
}


// Divide
template<typename T, typename U>
auto operator/(T a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(b.size());
  std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a/x; });
  return c;
}
template<typename T, typename U>
auto operator/(const vector<T>& a, U b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x/b; });
  return c;
}
template<typename T, typename U>
auto operator/(const vector<T>& a, const vector<U>& b) {
  vector<std::common_type_t<T,U>> c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::divides<T>());
  return c;
}

template<typename T, typename U>
auto operator/=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x/b; });
  return a;
}
template<typename T, typename U>
auto operator/=(vector<T>& a, const vector<U>& b) {
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::divides<T>());
  return a;
}


namespace vec {

  // Dot Product
  template<typename T, typename U>
  auto dotProduct(const vector<T>& a, const vector<U>& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  }


  // Cross Product
  template<typename T, typename U>
  auto crossProduct(const vector<T>& a, const vector<U>& b) {
    assert(a.size()==3);
    assert(b.size()==3);
    vector<std::common_type_t<T,U>> c(3);
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
  }


  // Sum
  template<typename T>
  auto sum(const vector<T>& a) {
    return std::accumulate(a.begin(), a.end(), T(0));
  }


  // Norm
  template<typename T>
  auto norm(const vector<T>& a) {
    return std::sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.0));
  }


  // Element-wise absolute value
  template<typename T>
  vector<T> abs(const vector<T>& a) {
    vector<T> b(a.size());
    std::transform(a.begin(), a.end(), b.begin(), [](T x){ return (x<0) ? -x : x; });
    return b;
  }

  // Element-wise square root
  template<typename T>
  vector<T> sqrt(const vector<T>& a) {
    vector<T> b(a.size());
    std::transform(a.begin(), a.end(), b.begin(), [](T x){ return std::sqrt(x); });
    return b;
  }

  // Element-wise exponent
  template<typename T, typename U>
  vector<T> pow(const vector<T>& a, U n) {
    vector<T> b(a.size());
    std::transform(a.begin(), a.end(), b.begin(), [n](T x){ return std::pow(x, n); });
    return b;
  }

  template<typename T, typename U>
  vector<T> pow(T a, const vector<U>& n) {
    vector<T> b(n.size());
    std::transform(n.begin(), n.end(), b.begin(), [a](U ni){ return std::pow(a, ni); });
    return b;
  }


  // Logical
  template<typename T>
  bool any(const std::vector<T>& a) {
    return std::any_of(a.begin(), a.end(), [](bool i){ return i; });
  }

  template<typename T>
  bool all(const std::vector<T>& a) {
    return std::all_of(a.begin(), a.end(), [](bool i){ return i; });
  }


  // Insert
  template<typename T>
  void insert_unique(std::vector<T>& vec, T value) {
    for (auto iter=vec.begin(); iter<vec.end(); iter++) {
      if (value <= *iter) {
        if (value != *iter) vec.insert(iter, value);
        return;
      }
    }
    vec.push_back(value);
  }

}

#endif

#ifndef VEC_HPP
#define VEC_HPP

#include <set>
#include <math.h>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <functional>


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
auto& operator+=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x+b; });
  return a;
}
template<typename T, typename U>
auto& operator+=(vector<T>& a, const vector<U>& b) {
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
auto& operator-=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x-b; });
  return a;
}
template<typename T, typename U>
auto& operator-=(vector<T>& a, const vector<U>& b) {
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
auto& operator*=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x*b; });
  return a;
}
template<typename T, typename U>
auto& operator*=(vector<T>& a, const vector<U>& b) {
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
auto& operator/=(vector<T>& a, U b) {
  std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x/b; });
  return a;
}
template<typename T, typename U>
auto& operator/=(vector<T>& a, const vector<U>& b) {
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


  // Product
  template<typename T>
  auto product(const vector<T>& a) {
    if (a.empty()) return T(0);
    return std::accumulate(a.begin(), a.end(), T(1), std::multiplies<T>());
  }


  // Norm
  template<typename T>
  auto norm(const vector<T>& a) {
    return std::sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.0));
  }


  // Root mean square
  template<typename T>
  auto rms(const vector<T>& a) {
    return std::sqrt(vec::dotProduct(a,a) / a.size());
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


  // Check if value is contained in vector
  template<typename T>
  bool isIn(const std::vector<T>& vec, T value) {
    std::set<T> set(vec.begin(), vec.end());
    bool isInVector = (set.find(value) != set.end());
    return isInVector;
  }


  // Element-wise logical
  template<typename T>
  std::vector<char> lessThan(const std::vector<T>& v, T s) {
    std::vector<char> result(v.size());
    std::transform(v.begin(), v.end(), result.begin(), [s](T vi) { return vi < s; });
    return result;
  }

  template<typename T>
  std::vector<char> lessThan(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<char> result(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T a, T b) { return a < b; });
    return result;
  }

  template<typename T>
  std::vector<char> greaterThan(const std::vector<T>& v, T s) {
    std::vector<char> result(v.size());
    std::transform(v.begin(), v.end(), result.begin(), [s](T vi) { return vi > s; });
    return result;
  }

  template<typename T>
  std::vector<char> greaterThan(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<char> result(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T a, T b) { return a > b; });
    return result;
  }


  template<typename T>
  vector<T> slice(const vector<T>& in, const vector<int>& index) {
    vector<T> slice;
    slice.reserve(index.size());
    for (int i: index) slice.push_back(in[i]);
    return slice;
  }


  // Get the sorted vector
  template<typename T>
  vector<T> sort(const vector<T>& in, vector<int>* index) {
    if (!index) {
      vector<T> out = in;
      std::sort(out.begin(), out.end());
      return out;
    }
    *index = vector<int>(in.size());
    std::iota(index->begin(), index->end(), 0);
    std::stable_sort(index->begin(), index->end(), [&in](int i1, int i2) {return in[i1] < in[i2];});
    vector<T> out(index->size());
    for (size_t i=0; i<out.size(); i++) {
      out[i] = in[(*index)[i]];

    }
    return out;
  }

  // Get the sorted unique vector
  template<typename T>
  vector<T> unique(const vector<T>& in, vector<int>* index) {
    vector<T> out = sort(in, index);
    if (!index) {
      out.erase(std::unique(out.begin(), out.end()), out.end());
      return out;
    }
    T value = out[0];
    auto iIdx = index->begin()+1;
    for (auto iOut=out.begin()+1; iOut!=out.end();) {
      if (*iOut == value) {
        iOut = out.erase(iOut);
        iIdx = index->erase(iIdx);
      } else {
        value = *iOut;
        iOut++;
        iIdx++;
      }
    }
    return out;
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


  // Create a vector with an incremented value
  template<typename T, typename T1, typename T2>
  vector<T> arange(T1 start, T2 stop, T step, bool inclusive) {
    // Get number of values in array
    float nFloat = (float)(stop - start) / step;
    int n;
    if (fmod(nFloat,1)==0) {
      n = (int)nFloat; // Ensure that nFloat rounds correctly if is an integer with machine error
      if (inclusive) n++;
    } else {
      n = ceil(nFloat);
    }
    // Create array
    vector<T> vec(n);
    T value = start;
    for (int i=0; i<n; i++) {
      vec[i] = value;
      value += step;
    }
    return vec;
  }

}

#endif

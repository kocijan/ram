//
// Created by dario on 02. 06. 2020..
//

#ifndef RAM_UINT128_T_HPP
#define RAM_UINT128_T_HPP

#include <iostream>

struct uint128_t {
  std::uint64_t low;
  std::uint64_t high;

  // Empty constructor
  inline uint128_t() = default;
  // Promotion constructor
  inline uint128_t(std::uint64_t rhs) {
    high = 0;
    low = rhs;
  }

  inline uint128_t(std::uint64_t l, std::uint64_t h): low(l), high(h) {}

  // Promotion assignment operator
  inline uint128_t& operator=(std::uint64_t rhs) {
    high = 0;
    low = rhs;
    return *this;
  }

  // Conversion operator
  inline explicit operator std::uint64_t() const { return low; }

  inline explicit operator std::uint32_t() const { return static_cast<std::uint32_t>(low); }

  // Implicit boolean conversion
  inline explicit operator bool() const { return high || low; }

  inline uint128_t operator~() {
    uint128_t tmp(*this);
    tmp.low = ~tmp.low;
    tmp.high = ~tmp.high;
    return tmp;
  }

  // Increment operators
  inline uint128_t& operator++() {
    if (++low == 0) {
      ++high;
    }
    return *this;
  }
  inline uint128_t operator++(int) {
    uint128_t tmp(*this);
    operator++();
    return tmp;
  }

  // Decrement operators
  inline uint128_t& operator--() {
    --low;
    if (low == ~0ull) {
      --high;
    }
    return *this;
  }
  inline uint128_t operator--(int) {
    uint128_t tmp(*this);
    operator--();
    return tmp;
  }

  // Compound Arithemetic operators
  inline uint128_t& operator+=(const uint128_t& rhs) {
    std::uint64_t low = this->low;
    this->low += rhs.low;
    this->high += rhs.high;
    if (this->low < low) {
      ++this->high;
    }
    return *this;
  }
  inline uint128_t& operator-=(const uint128_t& rhs) {
    std::uint64_t low = this->low;
    this->low -= rhs.low;
    this->high -= rhs.high;
    if (this->low > low) {
      --this->high;
    }
    return *this;
  }

  // Compound Bitwise operators
  inline uint128_t& operator&=(const uint128_t& rhs) {
    this->low &= rhs.low;
    this->high &= rhs.high;
    return *this;
  }
  inline uint128_t& operator|=(const uint128_t& rhs) {
    this->low |= rhs.low;
    this->high |= rhs.high;
    return *this;
  }
  inline uint128_t& operator^=(const uint128_t& rhs) {
    this->low ^= rhs.low;
    this->high ^= rhs.high;
    return *this;
  }
  inline uint128_t& operator<<=(unsigned int shift) {
    // Left-shifting a uint64_t by 64+ is undefined behavior
    if (shift < 64) {
      std::uint64_t high_low_bits = this->low >> (64 - shift);
      this->low <<= shift;
      this->high <<= shift;
      this->high |= high_low_bits;
    } else {
      this->high = this->low << (shift - 64);
      this->low = 0;
    }
    return *this;
  }
  inline uint128_t& operator>>=(unsigned int shift) {
    if (shift < 64) {
      std::uint64_t low_high_bits = this->high << (64 - shift);
      this->low >>= shift;
      this->high >>= shift;
      this->low |= low_high_bits;
    } else {
      this->low = this->high << (shift - 64);
      this->high = 0;
    }
    return *this;
  }
};

// Arithmetic operators
inline uint128_t operator+(uint128_t lhs, const uint128_t& rhs) {
  lhs += rhs;
  return lhs;
}
inline uint128_t operator-(uint128_t lhs, const uint128_t& rhs) {
  lhs -= rhs;
  return lhs;
}

// Bitwise operators
inline uint128_t operator&(uint128_t lhs, const uint128_t& rhs) {
  lhs &= rhs;
  return lhs;
}
inline uint128_t operator|(uint128_t lhs, const uint128_t& rhs) {
  lhs |= rhs;
  return lhs;
}
inline uint128_t operator^(uint128_t lhs, const uint128_t& rhs) {
  lhs ^= rhs;
  return lhs;
}
inline uint128_t operator<<(uint128_t lhs, unsigned int shift) {
  lhs <<= shift;
  return lhs;
}
inline uint128_t operator>>(uint128_t lhs, unsigned int shift) {
  lhs >>= shift;
  return lhs;
}

// Comparison operators
inline bool operator==(const uint128_t& lhs, const uint128_t& rhs) {
  return lhs.high == rhs.high && lhs.low == rhs.low;
}
inline bool operator!=(const uint128_t& lhs, const uint128_t& rhs) {
  return !operator==(lhs, rhs);
}
inline bool operator<(const uint128_t& lhs, const uint128_t& rhs) {
  return (lhs.high < rhs.high) || (lhs.high == rhs.high && lhs.low < rhs.low);
}
inline bool operator>(const uint128_t& lhs, const uint128_t& rhs) {
  return operator<(rhs, lhs);
}
inline bool operator<=(const uint128_t& lhs, const uint128_t& rhs) {
  return !operator>(lhs, rhs);
}
inline bool operator>=(const uint128_t& lhs, const uint128_t& rhs) {
  return !operator<(lhs, rhs);
}

// Compatibility with namespace std primitives
namespace std {
// Hash function (for std::hash)
template <>
struct hash<uint128_t> {
  // adapted from boost::hash_combine
  inline size_t operator()(const uint128_t& obj) const {
    size_t low = hash<int64_t>()(obj.low);
    size_t high =
        hash<int64_t>()(obj.high) + 0x9e3779b9 + (low << 6) + (low >> 2);
    return low ^ high;
  }
};

// C++ Stream operators
inline ostream& operator<<(ostream& lhs, const uint128_t& rhs) {
  char buffer[35];
  sprintf(buffer, "0x%016lx%016lx", rhs.high, rhs.low);
  lhs << buffer;
  return lhs;
}
}  // namespace std

#endif  // RAM_UINT128_T_HPP
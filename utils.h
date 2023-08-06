#pragma once
#include <climits>
#include <cstring>
#include <iostream>

#include "parlay/utilities.h"
using namespace std;

constexpr size_t BLOCK_SIZE = 1024;

template <typename ET>
inline bool atomic_compare_and_swap(ET *a, ET oldval, ET newval) {
  static_assert(sizeof(ET) <= 8, "Bad CAS length");
  if (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t *>(a), r_oval,
                                        r_nval);
  } else if (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t *>(a), r_oval,
                                        r_nval);
  } else {  // if (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t *>(a), r_oval,
                                        r_nval);
  }
}

template <class ET>
inline bool compare_and_swap(std::atomic<ET> *a, ET oldval, ET newval) {
  return a->load() == oldval &&
         atomic_compare_exchange_weak(a, &oldval, newval);
}

template <class ET>
inline bool compare_and_swap(ET *a, ET oldval, ET newval) {
  return (*a) == oldval && atomic_compare_and_swap(a, oldval, newval);
}

template <typename E, typename EV>
inline E fetch_and_add(E *a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
  return oldV;
}

template <typename E, typename EV>
inline void write_add(E *a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
}

template <typename ET, typename F>
inline bool write_min(ET *a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(ET *a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <class ET>
inline ET _hash(ET a) {
  if (sizeof(ET) == 4) {
    return parlay::hash32(a);
  } else if (sizeof(ET) == 8) {
    return parlay::hash64(a);
  } else {
    std::cout << "hash bad length: " << sizeof(ET) << std::endl;
    abort();
  }
}

template <class ET>
inline ET _hash_2(ET a) {
  if (sizeof(ET) == 4) {
    return parlay::hash32_2(a);
  } else if (sizeof(ET) == 8) {
    return parlay::hash64_2(a);
  } else {
    std::cout << "hash bad length: " << sizeof(ET) << std::endl;
    abort();
  }
}

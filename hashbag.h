#pragma once
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "utils.h"
using namespace std;
using namespace parlay;

constexpr size_t MIN_BAG_SIZE = 1 << 10;
constexpr int MAX_PROBES = 100;
constexpr int EXP_SAMPLES = 64;

template <class ET>
class hashbag {
 protected:
  size_t n;
  double load_factor;
  int pointer;
  int num_hash_bag;
  size_t pool_size;
  ET empty;
  sequence<size_t> bag_size;
  sequence<uint32_t> counter;
  sequence<ET> pool;
  void set_values(ET u, int &local_pointer, size_t &len, size_t &idx,
                  size_t &rate) {
    local_pointer = pointer;
    // size of hash bag
    len = bag_size[local_pointer] - bag_size[local_pointer - 1];
    // insert into a random position between [2^i, 2^(i+1)]
    idx = (_hash(u) & (len - 1)) + bag_size[local_pointer - 1];
    // reciprocal sample rate
    rate = max(1.0, floor(len * load_factor / EXP_SAMPLES));
  }
  void clear() {
    size_t len = bag_size[pointer];
    parallel_for(
        0, len,
        [&](size_t i) {
          if (pool[i] != empty) {
            pool[i] = empty;
          }
        },
        BLOCK_SIZE);
    for (int i = 0; i <= pointer; i++) {
      counter[i] = 0;
    }
    pointer = 1;
  }

 public:
  hashbag() = delete;
  hashbag(size_t _n, size_t _min_bag_size = MIN_BAG_SIZE,
          double _load_factor = 0.5, ET _empty = numeric_limits<ET>::max())
      : n(_n), load_factor(_load_factor), empty(_empty) {
    pointer = 1;
    size_t cur_size = _min_bag_size;
    pool_size = 0;
    bag_size = sequence<size_t>(1, 0);
    while (pool_size < n / load_factor) {
      bag_size.push_back(cur_size + bag_size.back());
      pool_size += cur_size;
      cur_size *= 2;
    }
    num_hash_bag = bag_size.size();
    counter = sequence<uint32_t>(num_hash_bag, 0);
    pool = sequence<ET>(pool_size, empty);
  }
  size_t get_bag_capacity() { return pool_size; }
  size_t insert(ET u) {
    int local_pointer;
    size_t len, idx, rate;
    set_values(u, local_pointer, len, idx, rate);
    while (idx % rate == 0) {
      bool succeed = false;
      uint32_t ret = counter[local_pointer];
      while (ret < EXP_SAMPLES && !succeed) {
        succeed = compare_and_swap(&counter[local_pointer], ret, ret + 1);
        ret = counter[local_pointer];
      }
      if (ret >= EXP_SAMPLES && local_pointer + 1 < num_hash_bag) {
        compare_and_swap(&pointer, local_pointer, local_pointer + 1);
      }
      if (local_pointer != pointer) {
        set_values(u, local_pointer, len, idx, rate);
      }
      if (succeed || local_pointer + 1 == num_hash_bag) {
        break;
      }
    }
    for (size_t i = 0;; i++) {
      if (compare_and_swap(&pool[idx], empty, u)) {
        break;
      }
      idx = (idx == bag_size[local_pointer] - 1 ? bag_size[local_pointer - 1]
                                                : idx + 1);
      if (i % MAX_PROBES == 0) {
        if (local_pointer != pointer) {
          set_values(u, local_pointer, len, idx, rate);
        }
      }
      if (i == len) {
        if (local_pointer == pointer) {
          compare_and_swap(&pointer, local_pointer, local_pointer + 1);
        }
        local_pointer = pointer;
        set_values(u, local_pointer, len, idx, rate);
      }
    }
    return idx;
  }

  template <typename Out_Seq>
  size_t pack_into(Out_Seq out) {
    size_t len = bag_size[pointer];
    auto pred =
        delayed_seq<bool>(len, [&](size_t i) { return pool[i] != empty; });
    size_t ret =
        pack_into_uninitialized(pool.cut(0, len), pred, make_slice(out));
    clear();
    return ret;
  }
};

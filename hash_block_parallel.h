#ifndef hash_block_parallel_hpp
#define hash_block_parallel_hpp

#include "utils.h"
using namespace std;
constexpr int BLOCK_SIZE = 32;

// auxiliary function for power x^p
int mypower(int x, int p) {
  int res = 1;
  for (int i = 0; i < p; i++) {
    res *= x;
  }
  return res;
}

// build a table
template <typename T>
void build(const T seq, parlay::sequence<parlay::sequence<int>> &table_seq) {
  size_t k = seq.size() / BLOCK_SIZE;
  // pre-computed power table [p^(BLOCK_SIZE), p^(2 * BLOCK_SIZE), ... p^(logk *
  // BLOCK_SIZE)]
  int LOG2_k = FASTLOG2(k);
  parlay::sequence<int> block_power_table;
  for (int i = 0; i < LOG2_k; i++) {
    block_power_table.push_back(mypower(PRIME_BASE, int(BLOCK_SIZE * (i + 1))));
  }
  table_seq.resize(LOG2_k + 1);
  // pre-computed log table [1, 2, 4, ..., logk]
  parlay::sequence<int> block_log_table;
  for (int i = 0; i < LOG2_k + 1; i++) {
    block_log_table.push_back(1 << i);
  }

  // for the first dim of the table
  parlay::sequence<int> aux_inside_block_table;
  for (int i = 0; i <= (int)(BLOCK_SIZE); i++) {
    aux_inside_block_table.push_back(mypower(PRIME_BASE, i));
  }
  table_seq[0].resize(k);
  parlay::parallel_for(0, k, [&](int j) {
    // table_seq[0][j] = hash_value(seq, j * BLOCK_SIZE, (j + 1) * BLOCK_SIZE -
    // 1);
    table_seq[0][j] = 0;
    for (int b_j = j * BLOCK_SIZE; b_j <= (int)((j + 1) * BLOCK_SIZE - 1);
         b_j++) {
      table_seq[0][j] +=
          (aux_inside_block_table[(j + 1) * BLOCK_SIZE - 1 - b_j] *
           (int)seq[b_j]);
    }
  });
  // for the remaining dims of the table
  for (int i = 1; i < LOG2_k + 1; i++) {
    table_seq[i].resize(k - block_log_table[i] + 1);
    parlay::parallel_for(0, k - block_log_table[i] + 1, [&](int j) {
      table_seq[i][j] = table_seq[i - 1][j] * block_power_table[i - 1] +
                        table_seq[i - 1][j + block_log_table[i - 1]];
    });
  }
}

int qpow(int n) {
  int ret = 1;
  int a = PRIME_BASE;
  while(n) {
    if(n&1) {
      ret = ret*a;
    }
    a = a*a;
    n >>= 1;
  }
  return ret;
}

// function to build the two tables. `n` is
// the length of the initial sequence size.
template <typename T>
size_t construct_table(T &A, T &B,
                       parlay::sequence<parlay::sequence<int>> &table_A,
                       parlay::sequence<parlay::sequence<int>> &table_B,
                       parlay::sequence<std::pair<int, int>> &pre_su_a,
                       parlay::sequence<std::pair<int, int>> &pre_su_b,
                       parlay::sequence<int> &auxiliary_single_power_table,
                       size_t n) {
  parlay::internal::timer t;
  // logn

  /**
   * standard block size
   */
  int BLOCK_SIZE_UPPER = FASTLOG2(n);

  /**
   * 32 / 64 ?
   */
  build(A, table_A);
  build(B, table_B);

  // build the suffix
  // auxiliary power table
  parlay::sequence<int> block_power_table;
  for (int i = 0; i < BLOCK_SIZE; i++) {
    block_power_table.push_back(mypower(PRIME_BASE, i));
  }

  size_t aux_size =
      std::max(table_A[0].size() * BLOCK_SIZE, table_B[0].size() * BLOCK_SIZE);
  //auxiliary_single_power_table.resize(aux_size);
  //auxiliary_single_power_table[0] = 1;
  //t.next("first");

  auxiliary_single_power_table = parlay::sequence<int>(aux_size, PRIME_BASE);
  auxiliary_single_power_table[0] = 1;

  scan_inplace(make_slice(auxiliary_single_power_table), parlay::multiplies<int>());

  //for (int i = 1; i < aux_size; i++) {
    //auxiliary_single_power_table[i] =
        //PRIME_BASE * auxiliary_single_power_table[i - 1];
  //};
  int a_actual_size = BLOCK_SIZE * int(A.size() / BLOCK_SIZE);
  int b_actual_size = BLOCK_SIZE * int(B.size() / BLOCK_SIZE);
  // prefix_a.resize(a_actual_size);
  // prefix_b.resize(b_actual_size);
  
  pre_su_a = parlay::sequence<std::pair<int, int>>::uninitialized(a_actual_size);
  pre_su_b = parlay::sequence<std::pair<int, int>>::uninitialized(b_actual_size);
  //pre_su_a.resize(a_actual_size);
  //pre_su_b.resize(b_actual_size);

  //t.next("table and resize");

  // for both prefix and suffix precomputing
  
  int num_blocks_a = a_actual_size / BLOCK_SIZE;
  int num_blocks_b = b_actual_size / BLOCK_SIZE;
  
  parlay::parallel_for(0, num_blocks_a, [&](size_t i) {
      int s = i * BLOCK_SIZE;
      int e = (i + 1) * BLOCK_SIZE;
      pre_su_a[s].first = A[s];
      for(int j = s + 1; j < e; j++) {
        pre_su_a[j].first = pre_su_a[j - 1].first * PRIME_BASE + A[j];
      }

      pre_su_a[e - 1].second = A[e - 1];
      for(int j = e - 2; j >= s; j--) {
        pre_su_a[j].second = pre_su_a[j + 1].second * PRIME_BASE + A[j];
      }
  });

  //t.next("second");

  parlay::parallel_for(0, num_blocks_b, [&](size_t i) {
      int s = i * BLOCK_SIZE;
      int e = (i + 1) * BLOCK_SIZE;
      pre_su_b[s].first = B[s];
      for(int j = s + 1; j < e; j++) {
        pre_su_b[j].first = pre_su_b[j - 1].first * PRIME_BASE + B[j];
      }

      pre_su_b[e - 1].second = B[e - 1];
      for(int j = e - 2; j >= s; j--) {
        pre_su_b[j].second = pre_su_b[j + 1].second * PRIME_BASE + B[j];
      }
  });

  //t.next("third");
  return BLOCK_SIZE;
}

template <typename T>
bool compare_lcp(int p, int q, int z,
                 const parlay::sequence<parlay::sequence<int>> &table_A,
                 const parlay::sequence<parlay::sequence<int>> &table_B,
                 const parlay::sequence<std::pair<int, int>> &S_A,
                 const parlay::sequence<std::pair<int, int>> &S_B,
                 const parlay::sequence<int> &aux_power_table, int t, const T &A,
                 const T &B) {
  // size_t t = S_A.size() / table_A[0].size(); // BLOCK_SIZE
  if (t == 0) {
    return false;
  }
  if ((p + (1 << z) * t + t) >= (int)(table_A[0].size() * t) ||
      (q + (1 << z) * t + t) >= (int)(table_B[0].size() * t)) {
    return false;
  }
  // if ((p + (1 << z) * t) >= (int)(A.size()) ||
  //     (q + (1 << z) * t) >= (int)(B.size())) {
  //   return false;
  // }
  int next_block_A = (p / t) + 1;
  int next_block_B = (q / t) + 1;
  size_t rest_A_size = next_block_A * t - p;
  size_t rest_B_size = next_block_B * t - q;
  int hash_a_v;
  int hash_b_v;

  if (p % t == 0) {
    hash_a_v =
        table_A[z][p / t] * aux_power_table[t] + table_A[0][p / t + (1 << z)];
  } else {
    hash_a_v = S_A[p].second * aux_power_table[(1 << z) * t + t] +
               (table_A[z][next_block_A]) * aux_power_table[t - rest_A_size] +
               S_A[p + (1 << z) * t + t].first;
  }

  if (q % t == 0) {
    hash_b_v =
        table_B[z][q / t] * aux_power_table[t] + table_B[0][p / t + (1 << z)];
  } else {
    hash_b_v = S_B[q].second * aux_power_table[(1 << z) * t + t] +
               (table_B[z][next_block_B]) * aux_power_table[t - rest_B_size] +
               S_B[q + (1 << z) * t + t].first;
  }
  if (hash_a_v == hash_b_v) return true;
  return false;
}

// function for query the lcp from A[p] and B[q]
template <typename T>
int block_query_lcp(int p, int q, const T &A, const T &B,
                    const parlay::sequence<parlay::sequence<int>> &table_A,
                    const parlay::sequence<parlay::sequence<int>> &table_B,
                    const parlay::sequence<std::pair<int, int>> &S_A,
                    const parlay::sequence<std::pair<int, int>> &S_B,
                    const parlay::sequence<int> &aux_power_table, int t) {
  // find the possible block range point (omit the offset first)
  if (p >= (int)(A.size()) || q >= (int)(B.size())) {
    return 0;
  }

  if (A[p] != B[q]) {
    return 0;
  }
  int x = 0;
  // size_t t = S_A.size() / table_A[0].size(); // BLOCK_SIZE
  while ((p + t * mypower(2, x) + t < (int)(table_A[0].size() * t)) &&
         (q + t * mypower(2, x) + t < (int)(table_B[0].size() * t))) {
    if (!compare_lcp(p, q, x, table_A, table_B, S_A, S_B, aux_power_table, t, A,
                     B)) {
      break;
    }
    x++;
  }
  int pp = p;
  int qq = q;

  if (x > 0) {
    pp = p + t * mypower(2, x - 1) + t;
    qq = q + t * mypower(2, x - 1) + t;
    int y = x - 1;
    while (y >= 0) {
      if (compare_lcp(pp, qq, y, table_A, table_B, S_A, S_B, aux_power_table, t,
                      A, B)) {
        pp += t * mypower(2, y) + t;
        qq += t * mypower(2, y) + t;
      }
      y--;
    }
  } else {
    pp = p;
    qq = q;
  }

  while (pp < (int)(A.size()) && qq < (int)(B.size()) &&
         (int(A[pp]) == int(B[qq]))) {
    pp++;
    qq++;
  }

  if (pp != p) {
    return (pp - p);
  } else {
    return 1;
  }
}

#endif

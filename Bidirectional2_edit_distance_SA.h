#include "edit_distance_parallel.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>
#include <chrono>
#include <queue>
#include <mutex>
#include <cmath>

#include "dc3.h"
#include "parlay/sequence.h"
#include "range_min.h"
#include "sparse_table_sequential.h"
#include "suffix_array_parallel.h"
#include "utils.h"

// #define COMPUTE_AVERAGE_LCP

size_t bidirectional2_edit_distance_SA(const parlay::sequence<uint32_t> &a,
                      const parlay::sequence<uint32_t> &b, double *building_tm,
                      bool use_DC3)
{
  Timer tmr;
#ifdef COMPUTE_AVERAGE_LCP
  std::atomic_uint64_t lcp_total = 0;
  std::atomic_uint64_t lcp_cnt = 0;
#endif
  int n = a.size(), m = b.size();
  auto c = parlay::sequence<uint32_t>(n + m);
  parlay::parallel_for(0, n, [&](int i)
                       { c[i] = a[i]; });
  parlay::parallel_for(0, m, [&](int i)
                       { c[i + n] = b[i]; });
  auto rank = parlay::sequence<unsigned int>();
  auto sa = parlay::sequence<unsigned int>();
  auto lcp = parlay::sequence<unsigned int>();
  if (use_DC3)
  {
    std::tie(rank, sa, lcp) = DC3(c);
  }
  else
  {
    std::tie(rank, sa, lcp) = suffix_array_large_alphabet(c);
  }
  auto rmq = range_min(lcp);
  auto GetLcp = [&](int i, int j) -> int
  {
    if (i == n || j == m)
      return 0;
    assert(0 <= i && i < n && 0 <= j && j < m);
    for (int k = 0; k < 8; k++)
    {
      if (i + k >= n || j + k >= m || a[i + k] != b[j + k])
      {
        return k;
      }
    }
    int l = rank[i], r = rank[j + n];
    if (l > r)
      std::swap(l, r);
    int id = rmq.query(l + 1, r);
    return std::min(lcp[id], (unsigned int)n - i);
  };
  *building_tm = tmr.elapsed();
  auto Diag = [&](int i, int j) { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  parlay::sequence<int> old_max(50, -1);

   //Jakob Engel Tomfoolery
  parlay::sequence<uint32_t> a_reversed = a;
  std::reverse(a_reversed.begin(), a_reversed.end());
  parlay::sequence<uint32_t> b_reversed = b;
  std::reverse(b_reversed.begin(), b_reversed.end());
  auto c_rev = parlay::sequence<uint32_t>(n + m);
  parlay::parallel_for(0, n, [&](int i)
                       { c_rev[i] = a_reversed[i]; });
  parlay::parallel_for(0, m, [&](int i)
                       { c_rev[i + n] = b_reversed[i]; });
  int n2 = a_reversed.size();
  int m2 = b_reversed.size();
  parlay::sequence<parlay::sequence<int>> Rev_table_s1;
  parlay::sequence<parlay::sequence<int>> Rev_table_s2;
  parlay::sequence<int> rev_logN1;
  parlay::sequence<int> rev_powerN1;
  //build_hash_table(a_reversed, b_reversed, Rev_table_s1, Rev_table_s2, rev_powerN1, rev_logN1);
  parlay::sequence<int> rev_max_row(n + m + 1, -1), rev_temp(n + m + 1);
  parlay::sequence<int> rev_old_max(50, -1);
  auto rank_rev = parlay::sequence<unsigned int>();
  auto sa_rev = parlay::sequence<unsigned int>();
  auto lcp_rev = parlay::sequence<unsigned int>();
  if (use_DC3)
  {
    std::tie(rank_rev, sa_rev, lcp_rev) = DC3(c_rev);
  }
  else
  {
    std::tie(rank_rev, sa_rev, lcp_rev) = suffix_array_large_alphabet(c_rev);
  }
  auto rmq_rev = range_min(lcp_rev);
  auto GetLcp_rev = [&](int i, int j) -> int
  {
    if (i == n || j == m)
      return 0;
    assert(0 <= i && i < n && 0 <= j && j < m);
    for (int k = 0; k < 8; k++)
    {
      if (i + k >= n || j + k >= m || a_reversed[i + k] != b_reversed[j + k])
      {
        return k;
      }
    }
    int l = rank_rev[i], r = rank_rev[j + n];
    if (l > r)
      std::swap(l, r);
    int id = rmq_rev.query(l + 1, r);
    return std::min(lcp_rev[id], (unsigned int)n - i);
  };

  auto get_index = [&](int index) -> int {
        return index < 0 ? old_max.size() + index : index;
    };

  auto get_index_rev = [&](int index) -> int {
        return index < 0 ? rev_old_max.size() + index : index;
    };


  // Above is the play area - Beware
  *building_tm = tmr.elapsed();
  max_row[Diag(0, 0)] = GetLcp(0, 0);
  // Jakob's værk
  rev_max_row[Diag(0, 0)] = GetLcp_rev(0, 0);

  // cout << "hashing: " << query_lcp(a, b, table_s1, table_s2, logN1, 1000,
  // 1000)
  //      << endl;
  // bfs for path
  // int k = 0;
  // This should probably be a parallel for loop in the next iteration! - Jakob
  // l and r should probably be the same variables and then just be called reversedly based on which direction
  // Maybe dont make two constructions, but make them overlap: Do this in the next iteration - Jakobs
  bool pathFound = false;
  int Total_length = n + m;

  int fk = 0;
  int bk = 0;

  for (;;)
  {
    if((bk + fk) > ((old_max.size()) * 0.75)){
        //size_t required_size = std::min(static_cast<size_t>(n + m + 1), static_cast<size_t>(old_max.size() * 3.5));
        size_t required_size = std::min(static_cast<size_t>(n + m + 1), (old_max.size() * 3));
        old_max.resize(required_size, -1);
        rev_old_max.resize(required_size, -1);
    }

    int i = 0;
    int ix = 2;
    parlay::parallel_for(i, ix, [&](int dir)
                         {
    if (dir == 0) {
      fk++;
      int l = Diag(0, std::min(fk, int(m)));
      int r = Diag(std::min(fk, int(n)), 0);
    //Forwards reached stop
      if (max_row[Diag(n, m)] == n) {return fk +bk;}  
    
    // find path
      for(int id = l; id < r+1; id++) {
        int t = -1;
        old_max[get_index(id-m)] = max_row[id];
        if (max_row[id] != -1) {
          int i = max_row[id];
          int j = i + m - id;
          if (i == n || j == m) {
            t = i;
          } else {
          int get_lcp =
              GetLcp(i + 1, j + 1);
          t = i + 1 + get_lcp;
          }
        }
        if (id > 0 && max_row[id - 1] != -1) {
          int i = max_row[id - 1];
          int j = i + m - id + 1;
          if (i == n) {
            t = n;
          } else {
            t = std::max(
              t, i + 1 + GetLcp(i + 1, j));
          }
        }
        if (id < n + m && max_row[id + 1] != -1) {
          int i = max_row[id + 1];
          int j = i + m - id - 1;
          if (j == m) {
            t = std::max(t, i);
          } else {
            t = std::max(
                t, i + GetLcp( i, j + 1));
          }
        }
        // assert(t <= n);
        temp[id] = t;
    }
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
  }
  else if (dir == 1){
    bk++;
    int l2 = Diag(0, std::min(bk, int(m2)));
    int r2 = Diag(std::min(bk, int(n2)), 0);

    //Backwards reached stop
      if (rev_max_row[Diag(n2, m2)] == n2) {return fk +bk;}
      
  //Here comes backwards boy! Dun dun dun! 
  // find path
    
    for (int id2 = l2; id2 < r2+1; id2++) {
      rev_old_max[get_index_rev(id2-m)] = rev_max_row[id2];
      int t2 = -1;
      if (rev_max_row[id2] != -1) {
        int i2 = rev_max_row[id2];
        int j2 = i2 + m2 - id2;
        if (i2 == n || j2 == m) {
          t2 = i2;
        } else {
          int get_lcp2 =
              GetLcp_rev(i2 + 1, j2 + 1);
          t2 = i2 + 1 + get_lcp2;
        }
      }
      if (id2 > 0 && rev_max_row[id2 - 1] != -1) {
        int i2 = rev_max_row[id2 - 1];
        int j2 = i2 + m2 - id2 + 1;
        if (i2 == n2) {
          t2 = n2;
        } else {
          t2 = std::max(
              t2, i2 + 1 + GetLcp_rev(i2 + 1, j2));
        }
      }
      if (id2 < n2 + m2 && rev_max_row[id2 + 1] != -1) {
        int i2 = rev_max_row[id2 + 1];
        int j2 = i2 + m2 - id2 - 1;
        if (j2 == m2) {
          t2 = std::max(t2, i2);
        } else {
          t2 = std::max(
              t2, i2 + GetLcp_rev( i2, j2 + 1));
        }
      }
      // assert(t <= n);
      rev_temp[id2] = t2;
    }
    parlay::parallel_for(l2, r2 + 1,
                         [&](int id) { rev_max_row[id] = std::min(rev_temp[id], id); });
    }
  });
  
  int l = Diag(0, std::min(fk, int(m)));
  int r = Diag(std::min(fk, int(n)), 0);
  bool found = false;
  parlay::parallel_for(l, r + 1, [&](int id){
    if (max_row[id] + rev_max_row[Total_length-id] >= n){
      if ((max_row[id] + rev_old_max[get_index_rev(Total_length-id-m)] >= n || old_max[get_index(id-m)] + rev_max_row[Total_length-id] >= n) && !found){
        found = true;
      }
      pathFound = true;}
    }); 
    if (found){
      bk--;
    }

    if (pathFound)
    {
      return fk + bk;
    }
    
  }
  return fk + bk;
}

  //Above is the play area - Beware


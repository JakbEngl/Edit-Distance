#include "hash_block_parallel.h"
#include <chrono>
#include <iostream>
#include <queue>
#include <mutex>
#include <cmath>
#include <limits>

using hash_T = int64_t;
// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int EditDistanceBlockHashBidirectional(const Seq &a, const Seq &b,
                                  double *building_tm) {
  parlay::internal::timer tt;
  // build sparse table

  /**
   * rocord of building time
   */
  Timer tmr;

  int n = a.size();
  int m = b.size();
  if (n == 0) return m;
  if (m == 0) return n;
  //parlay::sequence<parlay::sequence<hash_T>> table_A;
  //parlay::sequence<parlay::sequence<hash_T>> table_B;
  //parlay::sequence<parlay::sequence<HashString>> table_A;
  //parlay::sequence<parlay::sequence<HashString>> table_B;
  //parlay::sequence<std::pair<hash_T, hash_T>> S_A;
  //parlay::sequence<std::pair<hash_T, hash_T>> S_B;
  parlay::sequence<hash_r_b_T> table_A;
  parlay::sequence<hash_r_b_T> table_B;
  //construct_table(a, b, table_A, table_B, S_A, S_B, std::min(n, m));
  build_rolling_blk(a, b, table_A, table_B);
  auto Diag = [&](int i, int j) { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  parlay::sequence<int> rev_max_row(n + m + 1, -1), rev_temp(n + m + 1);
  //tt.next("building");
    //Jakob Engel Tomfoolery
  /*Seq a_reversed = a;
  std::reverse(a_reversed.begin(), a_reversed.end());
  Seq b_reversed = b;
  std::reverse(b_reversed.begin(), b_reversed.end());
  int n2 = a_reversed.size();
  int m2 = b_reversed.size();
  parlay::sequence<parlay::sequence<hash_T>> table_A_rev;
  parlay::sequence<parlay::sequence<hash_T>> table_B_rev;
  parlay::sequence<std::pair<hash_T, hash_T>> S_A_rev;
  parlay::sequence<std::pair<hash_T, hash_T>> S_B_rev;
  construct_table(a_reversed, b_reversed, table_A_rev, table_B_rev, S_A_rev, S_B_rev, std::min(n, m));*/
  *building_tm = tmr.elapsed();
  //max_row[Diag(0,0)] = block_query_lcp(0, 0, a, b, table_A, table_B, S_A, S_B);
  max_row[Diag(0,0)] = query_rolling_blk(a,b,table_A, table_B, 0,0);
  // Jakob's værk
  //rev_max_row[Diag(n, m)] = block_query_lcp(0, 0, a_reversed, b_reversed, table_A_rev, table_B_rev, S_A_rev, S_B_rev);
  rev_max_row[Diag(n, m)] = n - query_rolling_blk_rev(a,b,table_A, table_B, n-1, m-1);
  // cout << "hashing: " << query_lcp(a, b, table_s1, table_s2, logN1, 1000,
  // 1000)
  //      << endl;
  // bfs for path
  int k = 0;
  bool pathFound = false;
  bool exceed = false;
  //int Total_length = n + m;
  for (;;) {
    //Forwards reached stop
    if (max_row[Diag(n, m)] == n) break;
    //Backwards reached stop
    //if (rev_max_row[Diag(n2, m2)] == n2) break;  
    // find path
    k++;
    int l = Diag(0, std::min(k, int(m)));
    int r = Diag(std::min(k, int(n)), 0);
    int l2 = Diag(std::max(n-k, 0), m);
    int r2 = Diag(n, std::max(m-k,0));
    parlay::parallel_for(l, r + 1, [&](int id) {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          int get_lcp = query_rolling_blk(a,b,table_A, table_B, i+1, j+1);
              //block_query_lcp(i + 1, j + 1, a, b, table_A, table_B, S_A, S_B);
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
              t, i + 1 + query_rolling_blk(a,b,table_A, table_B, i+1, j));
              //block_query_lcp(i + 1, j, a, b, table_A, table_B, S_A, S_B));
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          t = std::max(
              t, i + query_rolling_blk(a,b,table_A, table_B, i, j+1));
              //block_query_lcp(i, j + 1, a, b, table_A, table_B, S_A, S_B));
        }
      }
      // assert(t <= n);
      temp[id] = t;
    });
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
    
  //Check in the forwards one
  if(r >= l2){
    parlay::parallel_for(l2, r + 1, [&](int id){
      if (max_row[id] >= rev_max_row[id] && rev_max_row[id] != -1 ){
        if (max_row[id] > rev_max_row[id] && rev_max_row[id] != -1){exceed = true;}
        pathFound = true;}
      });
  }
  if (pathFound){
      k = pathFound ? k*2 -1: k;
      if (exceed){
        k--;
      }
      break;}

  //Here comes backwards boy! Dun dun dun! 
  // find path

    parlay::parallel_for(l2, r2 + 1, [&](int id2) {
      //int t2 = -1;
      int t2 = std::numeric_limits<int>::max();
      if (rev_max_row[id2] != -1) {
        int i2 = rev_max_row[id2];
        int j2 = i2 + m - id2;
        if (i2 == 0 || j2 == 0) {
          t2 = i2;
        } else {
          int get_lcp2 = query_rolling_blk_rev(a,b,table_A, table_B, i2-1, j2-1);
              //block_query_lcp(i2+1, j2 + 1, a_reversed, b_reversed, table_A_rev, table_B_rev, S_A_rev, S_B_rev);
          t2 = i2 - 1 - get_lcp2;
        }
      }
      if (id2 > 0 && rev_max_row[id2 - 1] != -1) {
        int i2 = rev_max_row[id2 - 1];
        int j2 = i2 + m - id2 + 1;
        if (i2 == 0) {
          t2 = 0;
        } else {
          t2 = std::min(
              t2, i2 - query_rolling_blk_rev(a,b, table_A, table_B, i2, j2-1));
              //block_query_lcp(i2+1, j2, a_reversed, b_reversed, table_A_rev, table_B_rev, S_A_rev, S_B_rev));
        }
      }
      if (id2 < n + m && rev_max_row[id2 + 1] != -1) {
        int i2 = rev_max_row[id2 + 1];
        int j2 = i2 + m - id2 - 1;
        if (j2 == 0) {
          t2 = std::min(t2, i2);
        } else {
          t2 = std::min(
              t2, i2 - 1 - query_rolling_blk_rev(a,b, table_A, table_B, i2-1, j2));
              //block_query_lcp(i2, j2 + 1, a_reversed, b_reversed, table_A_rev, table_B_rev, S_A_rev, S_B_rev));
        }
      }
      // assert(t <= n);
      rev_temp[id2] = t2;
    });
    parlay::parallel_for(l2, r2 + 1,
                         [&](int id) { rev_max_row[id] = std::min(rev_temp[id], id); });
    

    if (r >= l2){
    //std::cout << "We are going wide" << std::endl;
    parlay::parallel_for(l2, r +1, [&](int id){
      //if (max_row[id] + rev_max_row[id] >= n){pathFound = true;}
      if (max_row[id] >= rev_max_row[id] && rev_max_row[id] != -1){
        if (max_row[id] > rev_max_row[id] && rev_max_row[id] != -1){exceed = true;}
        pathFound = true;
        }
      });
  }
    if (pathFound){
      k = pathFound ? k*2 : k;
      if (exceed){
        k --;
      }
      break;}
  }
  return k;
}
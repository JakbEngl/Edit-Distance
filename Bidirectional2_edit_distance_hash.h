#include "hash_lcp_parallel.h"
#include "rolling_hashing.h"
#include <limits>
// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int bidirectional2_edit_distance_hash(const Seq &a, const Seq &b, double *building_tm)
{
  // build sparse table

  /**
   * Record of building time
   */
  Timer tmr;
  int n = a.size();
  int m = b.size();
  if (n == 0)
  {
    return m;
  }
  if (m == 0)
  {
    return n;
  }
  parlay::sequence<int> logN1;
  parlay::sequence<int> powerN1;
  //parlay::sequence<parlay::sequence<int>> table_s1;
  //parlay::sequence<parlay::sequence<int>> table_s2;
  parlay::sequence<hash_r_T> table_s1;
  parlay::sequence<hash_r_T> table_s2;

  //build_hash_table(a, b, table_s1, table_s2, powerN1, logN1);
  build_rolling(a, b, table_s1, table_s2);
  auto Diag = [&](int i, int j)
  { return i - j + m; };
  /*auto Diag2 = [&](int i, int j){
    return i - j;
  };*/

  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  //parlay::sequence<int> old_max(n - m + 20, -1);
  parlay::sequence<int> old_max(50, -1);
 
  // Jakob Engel Tomfoolery
  /*Seq a_reversed = a;
  std::reverse(a_reversed.begin(), a_reversed.end());
  Seq b_reversed = b;
  std::reverse(b_reversed.begin(), b_reversed.end());
  int n2 = a_reversed.size();
  int m2 = b_reversed.size();
  parlay::sequence<parlay::sequence<int>> Rev_table_s1;
  parlay::sequence<parlay::sequence<int>> Rev_table_s2;
  parlay::sequence<int> rev_logN1;
  parlay::sequence<int> rev_powerN1;
  build_hash_table(a_reversed, b_reversed, Rev_table_s1, Rev_table_s2, rev_powerN1, rev_logN1);*/

  parlay::sequence<int> rev_max_row(n + m + 1, -1), rev_temp(n + m + 1);
  //parlay::sequence<int> rev_old_max(n - m + 20, -1);
  parlay::sequence<int> rev_old_max(50, -1);
  auto get_index = [&](int index) -> int {
        return index < 0 ? old_max.size() + index : index;
    };
 
  auto get_index_rev = [&](int index) -> int {
        return index < 0 ? rev_old_max.size() + index : index;
    };


  // Above is the play area - Beware
  *building_tm = tmr.elapsed();
  //max_row[Diag(0, 0)] = query_lcp(a, b, table_s1, table_s2, logN1, 0, 0);
  max_row[Diag(0, 0)] = query_rolling(a,b,table_s1, table_s2, 0, 0);
  // Jakob's værk
  //rev_max_row[Diag(n, m)] = query_lcp(a_reversed, b_reversed, Rev_table_s1, Rev_table_s2, rev_logN1, 0, 0);
  rev_max_row[Diag(n, m)] = n - query_rolling_rev(a,b,table_s1, table_s2, n-1, m-1);

  // cout << "hashing: " << query_lcp(a, b, table_s1, table_s2, logN1, 1000,
  // 1000)
  //      << endl;
  // bfs for path
  // int k = 0;
  // This should probably be a parallel for loop in the next iteration! - Jakob
  // l and r should probably be the same variables and then just be called reversedly based on which direction
  // Maybe dont make two constructions, but make them overlap: Do this in the next iteration - Jakobs
  bool pathFound = false;
  //int Total_length = n + m;

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
    //if (max_row[Diag(n, m)] == n) {return fk +bk;}  
    
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
          int get_lcp = query_rolling(a,b, table_s1, table_s2, i+1, j+1);
              //query_lcp(a, b, table_s1, table_s2, logN1, i + 1, j + 1);
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
              t, i + 1 + query_rolling(a,b,table_s1, table_s2, i+1, j));
              //query_lcp(a, b, table_s1, table_s2, logN1, i + 1, j));
          }
        }
        if (id < n + m && max_row[id + 1] != -1) {
          int i = max_row[id + 1];
          int j = i + m - id - 1;
          if (j == m) {
            t = std::max(t, i);
          } else {
            t = std::max(
                t, i + query_rolling(a,b,table_s1, table_s2, i, j+1));
                //query_lcp(a, b, table_s1, table_s2, logN1, i, j + 1));
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
    //int l2 = Diag(0, std::min(bk, int(m2)));
    //int r2 = Diag(std::min(bk, int(n2)), 0);
    int l2 = Diag(std::max(n-bk, 0), m);
    int r2 = Diag(n, std::max(m-bk,0));

    //Backwards reached stop
    //if (rev_max_row[Diag(n2, m2)] == n2) {return fk +bk;}
      
  //Here comes backwards boy! Dun dun dun! 
  // find path
    
    for (int id2 = l2; id2 < r2+1; id2++) {
      rev_old_max[get_index_rev(id2-m)] = rev_max_row[id2];
      //int t2 = -1;
      int t2 = std::numeric_limits<int>::max();
      if (rev_max_row[id2] != -1) {
        int i2 = rev_max_row[id2];
        int j2 = i2 + m - id2;
        if (i2 == 0 || j2 == 0) {
          t2 = i2;
        } else {
          int get_lcp2 = query_rolling_rev(a,b, table_s1, table_s2, i2-1, j2-1);
              //query_lcp(a_reversed, b_reversed, Rev_table_s1, Rev_table_s2, rev_logN1, i2 + 1, j2 + 1);
          t2 = i2 - 1 - get_lcp2;
        }
      }
      if (id2 > 0 && rev_max_row[id2 - 1] != -1) {
        int i2 = rev_max_row[id2 - 1];
        int j2 = i2 + m - id2 + 1;
        if (i2 ==0) {
          t2 = 0;
        } else {
          t2 = std::min(
              t2, i2 - query_rolling_rev(a,b,table_s1, table_s2, i2, j2-1));
              //query_lcp(a_reversed, b_reversed, Rev_table_s1, Rev_table_s2, rev_logN1, i2 + 1, j2));
        }
      }
      if (id2 < n + m && rev_max_row[id2 + 1] != -1) {
        int i2 = rev_max_row[id2 + 1];
        int j2 = i2 + m - id2 - 1;
        if (j2 == 0) {
          t2 = std::min(t2, i2);
        } else {
          t2 = std::min(
              t2, i2 - 1 - query_rolling_rev(a,b,table_s1, table_s2, i2-1, j2));
              //query_lcp(a_reversed, b_reversed, Rev_table_s1, Rev_table_s2, rev_logN1, i2, j2 + 1));
        }
      }
      // assert(t <= n);
      rev_temp[id2] = t2;
    }
    parlay::parallel_for(l2, r2 + 1,
                         [&](int id) { rev_max_row[id] = std::min(rev_temp[id], id); });
    }
  });
  //int l = Diag(0, std::min(fk, int(m)));
  int r = Diag(std::min(fk, int(n)), 0);
  int l2 = Diag(std::max(n-bk, 0), m);
  bool found = false;
  if (r >= l2){
    parlay::parallel_for(l2, r + 1, [&](int id){
      if (max_row[id] >= rev_max_row[id] && rev_max_row[id] != -1){
        if (((max_row[id] >= rev_old_max[get_index_rev(id-m)] && rev_old_max[get_index_rev(id-m)] != -1 ) || (old_max[get_index(id-m)] >= rev_max_row[id] && rev_max_row[id] != -1)) && !found){
          found = true;
        }
        pathFound = true;}
      });}
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
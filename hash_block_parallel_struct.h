/*#ifndef hash_block_parallel_struct_hpp
#define hash_block_parallel_struct_hpp

#include "utils.h"

using namespace std;

using hash_T = int64_t;
constexpr int BLOCK_SIZE = 32;
constexpr int GRANULARITY = 512;

struct HashString {
    hash_T hash_value;
    string str;

    HashString(hash_T hash, const string& s) : hash_value(hash), str(s) {}

    HashString(const HashString& other) : hash_value(other.hash_value), str(other.str) {}

    HashString(HashString&& other) noexcept : hash_value(other.hash_value), str(std::move(other.str)) {}

    HashString& operator=(const HashString& other) {
        if (this != &other) {
            hash_value = other.hash_value;
            str = other.str;
        }
        return *this;
    }

    HashString& operator=(HashString&& other) noexcept {
        if (this != &other) {
            hash_value = other.hash_value;
            str = std::move(other.str);
        }
        return *this;
    }
};



// auxiliary function for power x^p
// int mypower(int x, int p) {
// int res = 1;
// for (int i = 0; i < p; i++) {
// res *= x;
//}
// return res;
//}

int mypower(int a, int n) {
  int ret = 1;
  while (n) {
    if (n & 1) {
      ret = ret * a;
    }
    a = a * a;
    n >>= 1;
  }
  return ret;
}

// build a table
template <typename T>
void build(const T &seq,
           parlay::sequence<parlay::sequence<HashString>> &table_seq) {
    std::cout << "Building table..." << std::endl;
    size_t k = seq.size() / BLOCK_SIZE;
    int LOG2_k = FASTLOG2(k);
    
    std::cout << "k: " << k << std::endl;

    // Initialize the table to store HashString instances
    table_seq = parlay::tabulate(LOG2_k + 1, [&](size_t i) {
        return parlay::sequence<HashString>::uninitialized(k);
    });

    std::cout << "Table initialized." << std::endl;

    // Calculate hashes and store them in the table
    parlay::parallel_for(0, k, [&](int j) {
        hash_T hash_value = 0;
        hash_T p = 1;
        string str = "";
        for (int b_j = (j + 1) * BLOCK_SIZE - 1; b_j >= j * BLOCK_SIZE; b_j--) {
            hash_value += seq[b_j] * p; // Adjust hash calculation according to your needs
            str += seq[b_j]; // Construct the string
            p *= PRIME_BASE;
        }
        // Store the hash and the string in HashString instance
        table_seq[0][j] = HashString(hash_value, str);
    });
     std::cout << "Hashes calculated" << std::endl;

    // Calculate remaining dimensions of the table
    for (int i = 1; i < LOG2_k + 1; i++) {
        parlay::parallel_for(
            0, k - (1 << i) + 1,
            [&](int j) {
                // Calculate hash for the combined string
                hash_T hash_value = table_seq[i - 1][j].hash_value *
                                        mypower(PRIME_BASE, BLOCK_SIZE * (1 << (i - 1))) +
                                    table_seq[i - 1][j + (1 << (i - 1))].hash_value;
                // Combine strings
                string str = table_seq[i - 1][j].str + table_seq[i - 1][j + (1 << (i - 1))].str;
                // Store the hash and the string in HashString instance
                table_seq[i][j] = HashString(hash_value, str);
            },
            GRANULARITY
        );
    }
    std::cout << "Remaining dimensions calculated" << std::endl;

    // Optionally, you can print the table for debugging purposes
    // for (size_t i = 0; i < table_seq.size(); i++) {
    //     printf("table[%zu]: ", i);
    //     for (size_t j = 0; j < table_seq[i].size(); j++) {
    //         printf("(%lld, %s) ", table_seq[i][j].hash_value, table_seq[i][j].str.c_str());
    //     }
    //     printf("\n");
    // }
}


int qpow(int n) {
  int ret = 1;
  int a = PRIME_BASE;
  while (n) {
    if (n & 1) {
      ret = ret * a;
    }
    a = a * a;
    n >>= 1;
  }
  return ret;
}

// function to build the two tables. `n` is
// the length of the initial sequence size.
template <typename T>
void construct_table(T &A, T &B,
                     parlay::sequence<parlay::sequence<HashString>> &table_A,
                     parlay::sequence<parlay::sequence<HashString>> &table_B,
                     parlay::sequence<std::pair<hash_T, hash_T>> &pre_su_a,
                     parlay::sequence<std::pair<hash_T, hash_T>> &pre_su_b,
                     size_t n) {
    if (A.size() < BLOCK_SIZE || B.size() < BLOCK_SIZE) {
        return;
    }
    
    // Calculate the upper bound for block size
    int BLOCK_SIZE_UPPER = FASTLOG2(n);

    // Build tables for both sequences A and B
    std::cout << "Building table for sequence A." << std::endl;
    build(A, table_A);
    std::cout << "Building table for sequence B." << std::endl;
    build(B, table_B);

    // Initialize pre_su_a and pre_su_b to store HashString instances
    int a_actual_size = BLOCK_SIZE * int(A.size() / BLOCK_SIZE);
    int b_actual_size = BLOCK_SIZE * int(B.size() / BLOCK_SIZE);
    std::cout << "Initializing pre_su_a with size: " << a_actual_size << std::endl;
    std::cout << "Initializing pre_su_b with size: " << b_actual_size << std::endl;
    pre_su_a = parlay::sequence<std::pair<hash_T, hash_T>>::uninitialized(a_actual_size);
    pre_su_b = parlay::sequence<std::pair<hash_T, hash_T>>::uninitialized(b_actual_size);

    // Calculate prefix and suffix for sequence A
    int num_blocks_a = a_actual_size / BLOCK_SIZE;
    std::cout << "Calculating prefix and suffix for sequence A with " << num_blocks_a << " blocks." << std::endl;
    parlay::parallel_for(0, num_blocks_a, [&](size_t i) {
        int s = i * BLOCK_SIZE;
        int e = (i + 1) * BLOCK_SIZE;
        pre_su_a[s].first = A[s];
        pre_su_a[s].second = A[e - 1];
        for (int j = s + 1; j < e; j++) {
            pre_su_a[j].first = pre_su_a[j - 1].first * PRIME_BASE + A[j];
        }
        for (int j = e - 2; j >= s; j--) {
            pre_su_a[j].second = pre_su_a[j + 1].second + (A[j] * qpow(e - j - 1));
        }
    });
    std::cout << "Processing done for blocks of sequence A." << std::endl;

    // Calculate prefix and suffix for sequence B (The bad place)
    int num_blocks_b = b_actual_size / BLOCK_SIZE;
    std::cout << "Calculating prefix and suffix for sequence B with " << num_blocks_b << " blocks." << std::endl;
    parlay::parallel_for(0, num_blocks_b, [&](size_t i) {
        int s = i * BLOCK_SIZE;
        int e = (i + 1) * BLOCK_SIZE;
        pre_su_b[s].first = B[s];
        pre_su_b[s].second = B[e - 1];

        for (int j = s + 1; j < e; j++) {
            pre_su_b[j].first = pre_su_b[j - 1].first * PRIME_BASE + B[j];
        }
        for (int j = e - 2; j >= s; j--) {
            pre_su_b[j].second = pre_su_b[j + 1].second + (B[j] * qpow(e - j - 1));
        }
    });
    std::cout << "Processing done for blocks of sequence B." << std::endl;

    // Optionally, you can print the prefix and suffix for debugging purposes
    // for (size_t i = 0; i < pre_su_a.size(); i++) {
    //     printf("pre_su_a[%zu]: (%lld,%lld)\n", i, pre_su_a[i].first, pre_su_a[i].second);
    // }
    // for (size_t i = 0; i < pre_su_b.size(); i++) {
    //     printf("pre_su_b[%zu]: (%lld,%lld)\n", i, pre_su_b[i].first, pre_su_b[i].second);
    // }
}


template <typename T>
bool compare_lcp(int p, int q, int z,
                 const parlay::sequence<parlay::sequence<HashString>> &table_A,
                 const parlay::sequence<parlay::sequence<HashString>> &table_B,
                 const parlay::sequence<std::pair<hash_T, hash_T>> &S_A,
                 const parlay::sequence<std::pair<hash_T, hash_T>> &S_B,
                 const T &A, const T &B) {
    constexpr int t = BLOCK_SIZE;
    if (t == 0) {
        return false;
    }
    int size_A = A.size() / t * t;
    int size_B = B.size() / t * t;
    if ((p + (1 << z) * t + t) >= size_A || (q + (1 << z) * t + t) >= size_B) {
        return false;
    }

    hash_T hash_a_v;
    hash_T hash_b_v;

    int block_idx_A = p / t;
    int block_idx_B = q / t;
    int offset_A = p % t;
    int offset_B = q % t;

    // Calculate hash_a_v
    if (offset_A == 0) {
        hash_a_v = table_A[z][block_idx_A].hash_value * qpow(t) + table_A[0][block_idx_A + (1 << z)].hash_value;
    } else {
        hash_a_v = S_A[p].second * qpow((1 << z) * t + (t - offset_A)) +
                   table_A[z][block_idx_A + 1].hash_value * qpow(t - offset_A) +
                   S_A[p + (1 << z) * t + t - 1].first;
    }

    // Calculate hash_b_v
    if (offset_B == 0) {
        hash_b_v = table_B[z][block_idx_B].hash_value * qpow(t) + table_B[0][block_idx_B + (1 << z)].hash_value;
    } else {
        hash_b_v = S_B[q].second * qpow((1 << z) * t + (t - offset_B)) +
                   table_B[z][block_idx_B + 1].hash_value * qpow(t - offset_B) +
                   S_B[q + (1 << z) * t + t - 1].first;
    }
    return hash_a_v == hash_b_v;
}


// function for query the lcp from A[p] and B[q]
template <typename T>
int block_query_lcp(int p, int q, const T &A, const T &B,
                    const parlay::sequence<parlay::sequence<HashString>> &table_A,
                    const parlay::sequence<parlay::sequence<HashString>> &table_B,
                    const parlay::sequence<std::pair<hash_T, hash_T>> &S_A,
                    const parlay::sequence<std::pair<hash_T, hash_T>> &S_B) {
    constexpr int t = BLOCK_SIZE;
    if (A.size() < t || B.size() < t) {
        int pp = p;
        int qq = q;
        while (pp < (int)A.size() && qq < (int)B.size() && A[pp] == B[qq]) {
            pp++;
            qq++;
        }
        return pp - p;
    }
    // find the possible block range point (omit the offset first)
    if (p >= (int)(A.size()) || q >= (int)(B.size())) {
        return 0;
    }

    if (A[p] != B[q]) {
        return 0;
    }
    int x = 0;
    int size_A = A.size() / t * t;
    int size_B = B.size() / t * t;
    while ((p + (t << x) + t < size_A) && (q + (t << x) + t < size_B)) {
        if (!compare_lcp(p, q, x, table_A, table_B, S_A, S_B, A, B)) {
            break;
        }
        x++;
    }

    int pp = p;
    int qq = q;

    if (x > 0) {
        pp = p + (t << (x - 1)) + t;
        qq = q + (t << (x - 1)) + t;
        int y = x - 1;
        while (y >= 0) {
            if (compare_lcp(pp, qq, y, table_A, table_B, S_A, S_B, A, B)) {
                pp += (t << y) + t;
                qq += (t << y) + t;
            }
            y--;
        }
    } else {
        pp = p;
        qq = q;
    }

    int cnt = 0;
    while (pp < (int)(A.size()) && qq < (int)(B.size()) &&
           (A[pp] == B[qq])) {
        cnt++;
        pp++;
        qq++;
    }
    return pp - p;
}


#endif
*/
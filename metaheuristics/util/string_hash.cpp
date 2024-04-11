#pragma once

#include <chrono>
#include <algorithm>
#include <vector>
#include "Constants.cpp"

const int M1 = 1000015553, M2 = 1000028537;
thread_local std::uniform_int_distribution<int> d1(356, M1 - 1), d2(356, M2 - 1);

const long long B1 = d1(g), B2 = d2(g);
 
struct HashedString {
    int N;
    std::string S;
    std::vector<long long> H1, P1, H2, P2;
    explicit HashedString() : N(0) { }
    explicit HashedString(std::string const& _S) : N(size(_S)), S(_S), H1(N, _S[0]), P1(N, 1), H2(N, _S[0]), P2(N, 1) {
        for(int i = 1; i < N; i++) {
            P1[i] = (P1[i - 1] * B1) % M1, P2[i] = (P2[i - 1] * B2) % M2;
            H1[i] = (H1[i - 1] * B1 + S[i]) % M1, H2[i] = (H2[i - 1] * B2 + S[i]) % M2;
        }
    }
    long long get(int l, int r) {
        if(l == 0) return ((H1[r] << 30) ^ (H2[r]));
        long long R1 = (((H1[r] - (H1[l - 1] * P1[r - l + 1])) % M1) + M1) % M1;
        long long R2 = (((H2[r] - (H2[l - 1] * P2[r - l + 1])) % M2) + M2) % M2;
        return ((R1 << 30) ^ (R2));
    }
};
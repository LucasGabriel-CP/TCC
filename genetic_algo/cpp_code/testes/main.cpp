#include <bits/stdc++.h>

#define rand_i() uid_i(g)
thread_local std::mt19937 g(std::chrono::high_resolution_clock::now().time_since_epoch().count());
thread_local std::uniform_int_distribution<int> uid_i(0, std::numeric_limits<int>::max());             // range

int main() {
    int start = 3, n = 7;
    for (int i = start + n; ; i--) {
        std::cout << i << " " << i - n << '\n';
        if (i == start) break;
    }
    std::cout << '\n';
    for (int i = start - n; ; i++) {
        std::cout << i << " " << i + n << '\n';
        if (i == start) break;
    }
    std::cout << '\n';
}
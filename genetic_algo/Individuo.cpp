#pragma once

#include <iostream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "LKC.cpp"


struct Segtree {
	std::vector<long long> seg, lazy;
	int n, LOG;

    Segtree(int _n) {
        n = _n;
        seg.assign(n, 0ll);
        lazy.assign(n, 0ll);
        LOG = std::ceil(log2(n));
    }

	long long junta(long long a, long long b) {
		return a+b;
	}

	void push(int p, long long x, int tam, bool prop=1) {
		seg[p] += x*tam;
		if (prop and p < n) lazy[p] += x;
	}

	void up(int p) {
		for (int tam = 2; p /= 2; tam *= 2) {
			seg[p] = junta(seg[2*p], seg[2*p+1]);
			push(p, lazy[p], tam, 0);
		}
	}

	void prop(int p) {
		int tam = 1 << (LOG-1);
		for (int s = LOG; s; s--, tam /= 2) {
			int i = p >> s;
			if (lazy[i]) {
				push(2*i, lazy[i], tam);
				push(2*i+1, lazy[i], tam);
				lazy[i] = 0;
			}
		}
	}

	long long query(int a, int b) {
		long long ret = 0;
		for (prop(a+=n), prop(b+=n); a <= b; ++a/=2, --b/=2) {
			if (a%2 == 1) ret = junta(ret, seg[a]);
			if (b%2 == 0) ret = junta(ret, seg[b]);
		}
		return ret;
	}

	void update(int a, int b, int x) {
		int a2 = a += n, b2 = b += n, tam = 1;
		for (; a <= b; ++a/=2, --b/=2, tam *= 2) {
			if (a%2 == 1) push(a, x, tam);
			if (b%2 == 0) push(b, x, tam);
		}
		up(a2), up(b2);
	}
};


struct Individuo {
    int fitness=-1;
    std::vector<std::vector<int>> matrix; // L por T

    Individuo() {
        matrix.assign(L, std::vector<int>(T, -1));
    }

    void dorit(int t, std::vector<Segtree> &facility_seg) {
        int low = 0, high = V-1;
        std::vector<int> ans;
        while (low < high) {
            int mid = (low + high) / 2;
            std::vector<int> S;
            std::unordered_set<int> set_V;
            for (int i = 0; i < V; i++) set_V.insert(i);
            while (!set_V.empty()) {
                int x = *set_V.begin();
                S.push_back(x);
                for (int i = 0; i <= mid; i++) {
                    int v = graph[x][i].first;
                    if (set_V.find(v) != set_V.end()) {
                        set_V.erase(v);
                        for (int j = 0; j <= mid; j++) {
                            int z = graph[v][j].first;
                            if (set_V.find(z) != set_V.end()) {
                                set_V.erase(z);
                            }
                        }
                    }
                }
                if ((int)S.size() <= K) {
                    ans = S;
                    high = mid;
                }
                else {
                    low = mid + 1;
                }
            }
        }
        for (int i: ans) {
            facility_seg[i].update(t, t, 1);
        }
        mtx.lock();
        debug(ans);
        mtx.unlock();
    }


    int jorge_algos(
        int t, std::vector<int> &dp, std::vector<Segtree> &facility_seg,
        std::vector<std::pair<int, int>> &rec, std::vector<std::vector<bool>> &vis
    ) {
        if (t >= T) return 0;
        auto &d = dp[t];
        if (d != -1) return d;
        d = 0;
        for (int l = 0; l < L; l++) {
            int lim = std::min(T, t + center_types[l]);
            for (int v = 0; v < V; v++) {
                if (vis[t][v]) continue;
                int qnt = facility_seg[v].query(t, lim-1);
                int aux = jorge_algos(lim, dp, facility_seg, rec, vis) + qnt;
                if (aux > d)  {
                    rec[t] = {l, v};
                    d = aux;
                }
            }
        }
        return d;
    }


    void build() {
        std::vector<Segtree> facility_seg(V, Segtree(T));
        for (int t = 0; t < T; t++) {
            dorit(t, facility_seg);
        }

        std::vector<std::pair<int, int>> rec(T);
        std::vector<std::vector<bool>> vis(T, std::vector<bool>(V, false));
        for (int k = 0; k < K; k++) {
            std::vector<int> dp(T, -1);
            jorge_algos(0, dp, facility_seg, rec, vis);
            int at = 0;
            while(at < T) {
                auto [l, v] = rec[at];
                matrix[l][at] = v;
                int nxt = std::min(T, at + center_types[l]);
                facility_seg[v].update(at, nxt-1, -1);
                for (int t = at; t < nxt; t++) {
                    vis[t][v] = true;
                }
                at = std::min(T, at + center_types[l]);
            }
        }
        
        // assert(validate());
        get_fitness();
    }

    void fix_this() {
        for (int l = 0; l < L; l++) {
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    int lim = std::min(T, t + center_types[t]);
                    for (int i = t; i < lim; i++) {
                        matrix[l][t] = -1;
                    }
                }
            }
        }
    }

    void mutate() {
        if (rand() < 0.5) {
            int idx = rand_i() % L;
            int l = rand_i() % T, r = rand_i() % T;
            if (l > r) std::swap(l, r);

            for (; l < r; l++, r--) {
                std::swap(matrix[idx][l], matrix[idx][r]);
            }
        }
        else {
            int idx_1 = rand_i() % L, idx_2 = rand_i() % L;
            int l = rand_i() % T, r = rand_i() % T;
            if (l > r) std::swap(l, r);
            for (; l <= r; l++) {
                std::swap(matrix[idx_1][l], matrix[idx_2][l]);
            }
        }

        get_fitness();
    }

    static Individuo cross(Individuo const &p1, Individuo const &p2) {
        Individuo children_1, children_2;
        if (rand() < 0.5) {
            int lim = 0;
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < T; j++) {
                    if (j > lim) {
                        children_1.matrix[i][j] = p1.matrix[i][j];
                        children_2.matrix[i][j] = p2.matrix[i][j];
                    }
                    else {
                        children_1.matrix[i][j] = p2.matrix[i][j];
                        children_2.matrix[i][j] = p1.matrix[i][j];
                    }
                }
                lim++;
            }
        }
        else {
            for (int i = 0; i < L; i++) {
                if (rand() < 0.5) {
                    children_1.matrix[i] = p1.matrix[i];
                    children_2.matrix[i] = p2.matrix[i];
                }
                else {
                    children_1.matrix[i] = p2.matrix[i];
                    children_2.matrix[i] = p1.matrix[i];
                }
            }
        }

        // Select a random child to return
        if (rand() < 0.5) {
            children_1.fix_this();
            children_1.get_fitness();
            return children_1;
        }
        children_2.fix_this();
        children_2.get_fitness();
        return children_2;
    }

    bool validate() {        
        #ifdef PIZZA
            // Verify solution
        #endif

        return true;
    }

    int give_penalties() {
        int total_penalty = 0;

        return total_penalty;
    }

    void local_search() {
        std::vector<Segtree> open_facilities(L, Segtree(T)), used_nodes(V, Segtree(T));
        std::vector<int> open_time(T, 0);
        for (int l = 0; l < L; l++){
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    int lim = std::min(T, t + center_types[l]);
                    open_facilities[l].update(t, lim-1, 1);
                    used_nodes[matrix[l][t]].update(t, lim-1, 1);
                    open_time[t]++;
                }
            }
        }

        for (int t = 0; t < T; t++) {
            if (open_time[t] < K) {
                int l = rand_i() % L;
                int lhs = std::min(T, t + center_types[l]) - 1;
                int rhs = std::min(0, t - center_types[l] + 1);
                int qnt = open_facilities[l].query(lhs, rhs);
                if (!qnt) {
                    int v = rand_i() % V;
                    qnt = used_nodes[v].query(lhs, rhs);
                    if (!qnt) {
                        matrix[l][t] = v;
                        t = rhs;
                    }
                }
            }
        }
    }

    int get_fitness() {
        fitness = 0;
        std::vector<std::vector<int>> open_facilities(T);
        for (int l = 0; l < L; l++) {
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    for (int i = t; i < std::min(T, t + center_types[l]); i++) {
                        open_facilities[t].push_back(matrix[l][t]);
                    }
                }
            }
        }


        for (int t = 0; t < T; t++) {
            for (int client: clients[t]) {
                int mn = INF;
                for (int facility: open_facilities[t]) {
                    mn = std::min(mn, graph[client][facility].second);
                }
                fitness = std::max(mn, fitness);
            }
        }
        
        fitness += give_penalties();
        
        return fitness;
    }
    
    friend bool operator < (const Individuo lhs, const Individuo rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const Individuo &ind) {
        os << "Fitness: " << ind.fitness;
        return os;
    }
};

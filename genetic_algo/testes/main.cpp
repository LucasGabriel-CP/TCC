#include <bits/stdc++.h>

#define rand_i() uid_i(g)
thread_local std::mt19937 g(std::chrono::high_resolution_clock::now().time_since_epoch().count());
thread_local std::uniform_int_distribution<int> uid_i(0, std::numeric_limits<int>::max());             // range

int V, T, L, K, E;
std::vector<std::vector<std::pair<int, long long>>> graph, ordered_graph;      // dist(i, j)
std::vector<std::vector<int>> clients;          // Clientes no instante t
std::vector<int> center_types;                  // Tipos de centro

void read(std::string const &instance_path) {
    std::ifstream in;

    in.open(instance_path);

    in >> V >> T >> L >> K;

    graph.assign(V, std::vector<std::pair<int, long long>>(V));
    ordered_graph.assign(V, std::vector<std::pair<int, long long>>(V));
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            graph[i][j].first = j;
            in >> graph[i][j].second;
            ordered_graph[i][j] = graph[i][j];
        }
        std::sort(ordered_graph[i].begin(), ordered_graph[i].end(), [](auto a, auto b){
            return a.second < b.second;
        });
    }

    center_types.resize(L);
    std::cout << K << '\n';
    for (int i = 0; i < L; i++) {
        in >> center_types[i];
        std::cout << center_types[i] << ' ';
    }
    std::cout << '\n';

    in.close();
}

int main() {
    for (int i = 1; i <= 30; i++){
        std::string str = "/home/lucas/TCC/data/lck_instances/inst";
        str += std::to_string(i);
        str += ".txt";
        read(str);
    }
}
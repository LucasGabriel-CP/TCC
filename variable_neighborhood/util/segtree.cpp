#pragma once

struct SegTree{
    #define Mid ((l + r) >> 1)
    #define nil 0
    #define MAXT 155
    int n;
    int lazy[4*MAXT], seg[4*MAXT];
    int modify(int &a, int &b){
        return int(a + b);
    }
    
    void push(int id, int tam){
        if (lazy[id]){
            seg[id] += lazy[id] * (tam + 1);
            if (tam){
                lazy[2*id] += lazy[id];
                lazy[2*id+1] += lazy[id];
            }
            lazy[id] = 0;
        }
    }

    int query(int id, int l, int r, int x, int y){
        push(id, r-l);
        if(x > r || y < l)   return nil; //doens'int change ans
        if(l >= x && r <= y) return seg[id];
        int p1 = query(2*id, l, Mid, x, y);
        int p2 = query(2*id+1, Mid+1, r, x, y);
        return modify(p1, p2);
    }

    void update(int id, int l, int r, int x, int y, int val){
        push(id, r-l);
        if(x > r || y < l)  return;
        if(l >= x && r <= y){
            lazy[id] = val;
            push(id, r-l);
            return;
        }
        update(2*id, l, Mid, x, y, val);
        update(2*id+1, Mid+1, r, x, y, val);
        seg[id] = modify(seg[2*id], seg[2*id+1]);
    }

    SegTree() {
        n = MAXT;
        for (int i = 0; i < 4 * n; i++) {
            lazy[i] = 0;
            seg[i] = 0;
        }
    }

    SegTree(int _n){
        n = _n;
        for (int i = 0; i < 4 * n; i++) {
            lazy[i] = 0;
            seg[i] = 0;
        }
    }

    int query(int x, int y){ return query(1, 0, n-1, x, y); }
    void update(int x, int y, int val){ update(1, 0, n-1, x, y, val); }
};
#pragma once

#include <iostream>
#include <math.h>
#include "Constants.cpp"
#define pi acos(-1.0)

double DEG_to_RAD(double d){ return d*pi/180.0; }
double RAD_to_DEG(double r){ return r*180.0/pi; }

struct Point{
    int x = 0, y = 0;
    constexpr Point(){ }
    constexpr Point(int x_, int y_) : x(x_), y(y_){ }
    friend std::ostream &operator << (std::ostream &os, Point& p){
        os << p.x << ' ' << p.y;
        return os;
    }
    friend std::istream &operator >>(std::istream &is, Point &p){
        is >> p.x >> p.y;
        return is;
    }
    friend bool operator < (const Point& lhs, const Point& rhs){
        if (lhs.x != rhs.x) return lhs.x < rhs.x;
        return lhs.y < rhs.y;
    }
    friend bool operator > (const Point& lhs, const Point& rhs){
        if (lhs.x != rhs.x) return lhs.x > rhs.x;
        return lhs.y > rhs.y;
    }
    friend bool operator == (const Point& lhs, const Point& rhs){
        return (lhs.x == rhs.x) && (lhs.y == rhs.y);
    }
    friend long long operator ^ (const Point& lhs, const Point& rhs) {
        return 1ll * lhs.x * rhs.y - 1ll * lhs.y * rhs.x;
    }
    friend long long operator*(const Point& a, const Point& b)  {
        return 1ll*a.x * b.x + 1ll*a.y * b.y;
    }
    Point operator+(const Point &p) { return {x + p.x, y + p.y}; }
    Point operator-(const Point &p) { return {x - p.x, y - p.y}; }
    Point operator*(int d) { return {x * d, y * d}; }
    Point operator/(int d) { return {x / d, y / d}; }
};

double dist(Point const &a, Point const &b) {
    return hypot(a.x - b.x, a.y - b.y);
}


bool intersect(Point p, Point r, Point q, Point s) {
    Point pq = q - p;
    Point pr = r - p;
    Point qs = s - q;
    long long c1 = pr ^ qs;
    long long c2 = pq ^ pr;
    long long c3 = pq ^ qs;
    if (!c1) {
        if (c2) return 0;
        double t0 = (double)(pq * pr) / (double)(pr * pr);
        double t1 = t0 + (double)(qs * pr) / (double)(pr * pr);
        if (f_cmp(t0, t1) > 0) std::swap(t0, t1);
        return f_cmp(0, t0) < 0 && f_cmp(t1, 1) < 0;
    }
    double t = (double)c3 / (double)c1;
    double u = (double)c2 / (double)c1;
    return f_cmp(0, t) < 0 && f_cmp(t, 1) < 0 && f_cmp(0, u) < 0 && f_cmp(u, 1) < 0;
}
#pragma once

#include <iostream>
#include <math.h>
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
    friend int operator ^ (const Point& lhs, const Point& rhs) {
        return lhs.x * rhs.y - lhs.y * rhs.x;
    }
    Point operator+(const Point &p) { return {x + p.x, y + p.y}; }
    Point operator-(const Point &p) { return {x - p.x, y - p.y}; }
    Point operator*(int d) { return {x * d, y * d}; }
    Point operator/(int d) { return {x / d, y / d}; }
};

float dist(Point const &a, Point const &b) {
    return hypot(a.x - b.x, a.y - b.y);
}

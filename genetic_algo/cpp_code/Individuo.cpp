#pragma once

#include <iostream>
#include <random>
#include <algorithm>

#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"

struct Individuo {
    double fitness=-1;

    Individuo() {
        
    }

    void build() {
        
        // assert(validate());
        get_fitness();
    }

    void fix_this() {

    }

    void mutate() {

        get_fitness();
    }

    static Individuo cross(Individuo const &p1, Individuo const &p2) {
        Individuo children_1, children_2;
        if (rand() < 0.5) {
            // Single point
        }
        else {
            // Multiple points
        }

        // if (!children_1.validate()) {
        //     children_1.fix_this();
        // }
        // if (!children_2.validate()) {
        //     children_2.fix_this();
        // }

        // Select a random child to return
        if (rand() < 0.5) {
            children_1.get_fitness();
            return children_1;
        }
        children_2.get_fitness();
        return children_2;
    }

    bool validate() {        
        #ifdef PIZZA
            // Verify solution
        #endif

        return true;
    }

    double give_penalties() {

        double total_penalty = 0.0;
        // Get penalties from Constants.cpp

        return total_penalty;
    }

    double get_fitness() {
        fitness = 0.0;

        // Calculate fitness based on problem formulation

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

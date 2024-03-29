#pragma once

struct Cable {
    int id, capacity, max_usage;
    double price;
    Cable() { }
    Cable(int _id, int _capacity, double _price, int _max_usage) {
        id = _id;
        capacity = _capacity;
        price = _price;
        max_usage = _max_usage;
    }

    double get_cable() {
        if (max_usage) {
            max_usage--;
            return price;
        }
        return -1.;
    }
};

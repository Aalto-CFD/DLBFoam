#pragma once

#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

#include "DynamicList.H"
#include "chemistryLoad.H"

namespace Foam {

struct Node;
using node_ptr = std::shared_ptr<Node>;

struct Node {

    // chemistryLoad() = default;
    Node(const chemistryLoad& load)
        : value(load) {}

    chemistryLoad         value;
    node_ptr              parent;
    std::vector<node_ptr> children;
};

struct loadTree {

    // A utility function to create a new N-ary tree node
    static node_ptr new_node(const chemistryLoad& value) {
        node_ptr temp(new Node(value));
        temp->value = value;
        return temp;
    }

    static double compute_send_value(const node_ptr& parent_candidate, const chemistryLoad& load) {

        // This is too arbitrary
        return 0.3 * parent_candidate->value.value; // TODO!
    }

    static void add_child(std::vector<node_ptr>& candidates, const chemistryLoad& load) {

        node_ptr correct_parent;

        double largest = 0.0;
        for (const auto candidate : candidates) {

            if (candidate->value > largest) {
                largest        = candidate->value.value; // TODO FIX!
                correct_parent = candidate;
            }
        }

        double send_value = compute_send_value(correct_parent, load);
        auto   node       = new_node(load);

        correct_parent->value -= send_value;
        node->value += send_value;
        node->parent = correct_parent;
        correct_parent->children.push_back(node);
    }

    static std::pair<std::vector<chemistryLoad>, std::vector<chemistryLoad>>
    divide(const DynamicList<chemistryLoad>& in, double treshold) {

        std::vector<chemistryLoad> big, small;

        for (const auto& v : in) {
            if (v > treshold) {
                big.push_back(v);
            } else {
                small.push_back(v);
            }
        }

        std::sort(big.begin(), big.end());
        std::sort(small.begin(), small.end());

        return {big, small};
    }


    static void print_tree(const node_ptr& x, std::vector<bool> flag, int depth = 0, bool isLast = false) {
        if (x == NULL) return;

        for (int i = 1; i < depth; ++i) {

            if (flag[i] == true) {
                std::cout << "| "
                          << " "
                          << " "
                          << " ";
            }

            else {
                std::cout << " "
                          << " "
                          << " "
                          << " ";
            }
        }

        if (depth == 0)
            std::cout << x->value.value << '\n';

        else if (isLast) {
            std::cout << "+--- " << x->value.value << '\n';

            flag[depth] = false;
        } else {
            std::cout << "+--- " << x->value.value << '\n';
        }

        int it = 0;
        for (auto i = x->children.begin(); i != x->children.end(); ++i, ++it) {
            print_tree(*i, flag, depth + 1, it == (x->children.size()) - 1);
        }
        flag[depth] = true;
    }
};

} // namespace Foam
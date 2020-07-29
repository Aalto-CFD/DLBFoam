#pragma once

#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

#include "DynamicList.H"
#include "chemistryLoad.H"
#include "Pstream.H"

namespace Foam {

struct Node;
using node_ptr = std::shared_ptr<Node>;

struct Node {

    // chemistryLoad() = default;
    Node(const chemistryLoad& load)
        : balanced_load(load), original_load(load) {}

    chemistryLoad         balanced_load;
    chemistryLoad         original_load;
    
    node_ptr              parent;
    std::vector<node_ptr> children;
};

struct loadTree {

    // A utility function to create a new N-ary tree node
    static node_ptr new_node(const chemistryLoad& value) {
        node_ptr temp(new Node(value));
        temp->balanced_load = value;
        return temp;
    }

    static void add_child(const node_ptr& parent, const node_ptr& child) {

        child->parent = parent;
        parent->children.push_back(child);
    }


    static void set_balanced_to_original(node_ptr& node){

        if (node) {
            node->balanced_load = node->original_load;

            for (auto child : node->children){
                set_balanced_to_original(child);
            }

        }

    }


    static node_ptr find(node_ptr& node, int rank) {

        if (node == nullptr) {
            return nullptr;
        }

        if (node->balanced_load.rank == rank) {return node;}


        for (auto child : node->children){
            auto target = find(child, rank);
            if (target != nullptr) {
                return target;
            }
        }


        return nullptr;

    }

    template <class Finder, class SendValue>
    static void add_child(std::vector<node_ptr>& candidates,
                          const chemistryLoad&   load,
                          const Finder&          f,
                          const SendValue&       c) {

        node_ptr correct_parent = f(candidates);

        if (correct_parent) {

            scalar send_value = c(correct_parent, load); // compute_send_value(correct_parent,
                                                         // load);
            auto node = new_node(load);

            correct_parent->balanced_load -= send_value;
            node->balanced_load += send_value;

            add_child(correct_parent, node);
        }
    }

    /*

    static void add_child(std::vector<node_ptr>& candidates, const chemistryLoad& load) {

        node_ptr correct_parent;

        scalar largest = 0.0;
        for (const auto candidate : candidates) {

            if (candidate->value > largest) {
                largest        = candidate->value.value; // TODO FIX!
                correct_parent = candidate;
            }
        }

        scalar send_value = compute_send_value(correct_parent, load);
        auto   node       = new_node(load);

        correct_parent->value -= send_value;
        node->value += send_value;
        node->parent = correct_parent;
        correct_parent->children.push_back(node);
    }
    */


    static void print(const node_ptr& root, int n_nodes = Pstream::nProcs()){

        std::vector<bool> flag(n_nodes, true); 
        print_tree(root, flag);
    }

private:

    static void
    print_tree(const node_ptr& x, std::vector<bool> flag, int depth = 0, bool isLast = false) {
        if (x == NULL) return;

        for (int i = 1; i < depth; ++i) {

            if (flag[i] == true) {
                Info << "| "
                          << " "
                          << " "
                          << " ";
            }

            else {
                Info << " "
                          << " "
                          << " "
                          << " ";
            }
        }

        if (depth == 0)
            Info << x->balanced_load.value << '\n';

        else if (isLast) {
            Info << "+--- " << x->balanced_load.value << '\n';

            flag[depth] = false;
        } else {
            Info << "+--- " << x->balanced_load.value << '\n';
        }

        size_t it = 0;
        for (auto i = x->children.begin(); i != x->children.end(); ++i, ++it) {
            print_tree(*i, flag, depth + 1, it == (x->children.size()) - 1);
        }
        flag[depth] = true;
    }
};

} // namespace Foam
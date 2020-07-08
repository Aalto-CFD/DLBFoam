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

    static void add_child(const node_ptr& parent, const node_ptr& child) {

        child->parent = parent;
        parent->children.push_back(child);
    }

    static node_ptr inorder(const node_ptr& node, int rank) {


        if (node == nullptr) {
            return nullptr;
        }

        if (node->value.rank == rank) {return node;}


        //this doesnt work because the return value is ignored

//        for (auto i = node->children.begin(); i != node->children.end(); ++i, ++it) {
//            auto asd = inorder(*i, rank, it == (node->children.size()) - 1);
//        }







        /*

        // Total children count
        size_t total = node->children.size();

        std::cout << total << std::endl;



        // All the children except the last
        for (size_t i = 0; i < total - 1; ++i){
            inorder(node->children[i], rank);
        } 

        */

        //if (node->value.rank == rank) {
        //    return node;
        //}


        // Print the current node's data
        //std::cout << node->value.value << " ";

        // Last child
        //inorder(node->children[total - 1], rank);


    }

    template <class Finder, class SendValue>
    static void add_child(std::vector<node_ptr>& candidates,
                          const chemistryLoad&   load,
                          const Finder&          f,
                          const SendValue&       c) {

        node_ptr correct_parent = f(candidates);

        if (correct_parent) {

            double send_value = c(correct_parent, load); // compute_send_value(correct_parent,
                                                         // load);
            auto node = new_node(load);

            correct_parent->value -= send_value;
            node->value += send_value;

            add_child(correct_parent, node);
        }
    }

    /*

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
    */

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

    static void
    print_tree(const node_ptr& x, std::vector<bool> flag, int depth = 0, bool isLast = false) {
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

        size_t it = 0;
        for (auto i = x->children.begin(); i != x->children.end(); ++i, ++it) {
            print_tree(*i, flag, depth + 1, it == (x->children.size()) - 1);
        }
        flag[depth] = true;
    }
};

} // namespace Foam
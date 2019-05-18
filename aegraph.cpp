// Copyright 2019 Luca Istrate, Danut Matei
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        // Add the current edge to existing paths
        auto r = subgraphs[i].possible_double_cuts();
        for (auto& v : r) {
            v.insert(v.begin(), i);
        }

        // Check if we can "double cut" at the current son
        if (subgraphs[i].size() == 1 && subgraphs[i].num_atoms() == 0) {
            r.push_back({i});
        }
        paths.insert(paths.end(), r.begin(), r.end());
    }

    return paths;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    auto pos_remove = where[0];
    AEGraph copied(repr());

    // Replace the indicated subgraph with its new form
    if (where.size() > 1) {
        // The node needed to be removed is further down
        std::vector<int> new_where(where.begin() + 1, where.end());
        auto new_subgraph = copied.subgraphs[pos_remove].double_cut(new_where);
        copied.subgraphs[pos_remove] = new_subgraph;
    } else {
        // Delete the subgraph
        auto node = subgraphs[pos_remove].subgraphs[0];
        copied.subgraphs.erase(copied.subgraphs.begin() + pos_remove);

        // Copy the atoms and subgraphs from the deleted node
        for (const auto& sub : node.subgraphs) {
            copied.subgraphs.push_back(sub);
        }
        for (const auto& atom : node.atoms) {
            copied.atoms.push_back(atom);
        }
    }

    return copied;
}


std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();
    int len = size();

    //  Check if we can erase subgraphs
    for (int i = 0; i < len_subgraphs; ++i) {
        if (level % 2) {
	    	paths.push_back({i});
	    }

        auto r = subgraphs[i].possible_erasures(level + 1);
    	for (auto& v : r) {
            v.insert(v.begin(), i);
            if (subgraphs[i].size() > 1 || level == -1) {
        		paths.push_back(v);
            }
        }
    }

    // Check if we can erase atoms
    for (int i = len_subgraphs; i < len; ++i) {
    	if (level % 2) {
	    	paths.push_back({i});
	    }
    }

    return paths;
}


AEGraph AEGraph::erase(std::vector<int> where) const {
    AEGraph tmp(repr());

    // Remove operation
	if (where.size() == 1) {
		// Delete atom
		if (where[0] >= num_subgraphs()){
			tmp.atoms.erase(tmp.atoms.begin() + where[0] - num_subgraphs());
		} else {
			// Delete subgraph
			tmp.subgraphs.erase(tmp.subgraphs.begin() + where[0]);
		}
	} else {
		// Continue search until the end of sequence
		std::vector<int> nextWhere(where.begin() + 1, where.end());
		auto newSubgraph = tmp.subgraphs[where[0]].erase(nextWhere);
		tmp.subgraphs[where[0]] = newSubgraph;
	}

    return tmp;
}


std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();
    int sons = size();

    for (int i = 0; i < sons; i++) {
        // Get all the paths to a subgraph similar to the current son
        std::vector<std::vector<int>> r;
        if (i < len_subgraphs) {
            r = get_paths_to(subgraphs[i]);
        } else {
            r = get_paths_to(atoms[i - len_subgraphs]);
        }

        // Remove the direct path to the son
        r.erase(remove(r.begin(), r.end(), std::vector<int>{i}), r.end());

        // Add the new deiteration positions
        paths.insert(paths.end(), r.begin(), r.end());
    }

    return paths;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    return erase(where);
}


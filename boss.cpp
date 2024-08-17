#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <string>
#include <cstring>


bool in_range(int i, int min, int max) {
    if (i >= min && i <= max) return true;
    return false;
}

int min_in_range(int i, int j, int min, int max) {
    if (in_range(i, min, max) && in_range(j, min, max)) return std::min(i, j);
    else if (in_range(i, min, max) && !in_range(j, min, max)) return i;
    else if (!in_range(i, min, max) && in_range(j, min, max)) return j;
    else return -1;
}

int max_in_range(int i, int j, int min, int max) {
    if (in_range(i, min, max) && in_range(j, min, max)) return std::max(i, j);
    else if (in_range(i, min, max) && !in_range(j, min, max)) return i;
    else if (!in_range(i, min, max) && in_range(j, min, max)) return j;
    else return -1;
}

int to_int(char c) {

	switch(c){
		case '$': case '%': return 0; break;
		case 'A': case 'a': return 1; break;
		case 'C': case 'c': return 2; break;
		case 'G': case 'g': return 3; break;
		case 'T': case 't': return 4; break;
		default: break;

	}

	return 1;
}


class BOSS {
    /*
        Class to represent a de Bruijn graph
    */
public:
    int k;                                  // Lenght of node labels
    int m;                                  // Number of edges (m = |L| = |W|)
    int n;                                  // Number of nodes (n = L.count(1))

    sdsl::rrr_vector<> L;                   // Bit vector used to represent whether an edge is the last edge exiting a node
    sdsl::wt_hutu<sdsl::rrr_vector<63>> W;  // Stores the labels of the edges
    sdsl::int_vector<> F;                   // Stores the first occurences of the last character of the labels of the nodes

private:
    sdsl::rrr_vector<>::rank_0_type* m_L_rank_0;
    sdsl::rrr_vector<>::rank_1_type* m_L_rank_1;
    sdsl::rrr_vector<>::select_0_type* m_L_select_0;
    sdsl::rrr_vector<>::select_1_type* m_L_select_1;
    int m_p;

public:
    BOSS(const std::string& file_path, int k) : k(k) {        
        /*
            Constructs the de Bruijn graph loading the string T from disk

            file_path: path to input file
            k: lenght of k-mers
        */

        // Construct the Compressed Suffix Array
        sdsl::csa_wt<> csa;

        sdsl::construct(csa, file_path, 1);

        // std::cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << std::endl;
        // csXprintf(std::cout, "%2I %2S %3s %3P %2p %3B   %:3T", csa); 
        // std::cout << std::endl;

        // We have to substract the k terminator characters added at the beginning        
        this->m = csa.size() - k;
        std::cout << "starting construction with m = " << this->m << '\n';

        // Temporary storage
        sdsl::bit_vector L_bv(this->m, 1);
        std::string W_str = "";
        sdsl::int_vector<> F_iv(csa.sigma-1, 0);

        // Variables needed for the loop
        int pos;
        std::string last_label(k, ' ');
        sdsl::bit_vector chars_seen(csa.sigma-1, 0);
        std::string label;

        // std::cout << "starting for loop\n";
        for (int i = k; i < csa.size(); i++) {

            // Read node label
            pos = csa[i];
            label = sdsl::extract(csa, pos, pos+k-1);
            reverse(label.begin(), label.end());

            // Read edge
            char edge = csa.bwt[i];

            // Check needed for ending $_M
            if (edge == '\0') {
                edge = '$';
                this->m_p = i - k;
            }

            // L -> store 0 for the first edges of each node
            if (label == last_label) {
                L_bv[i-k-1] = 0;
            }

            // W -> we use lowercase chars for flagged characters
            if (label.substr(1, k-1) == last_label.substr(1, k-1)) {
                // Turn an edge into lowercase if it needs to be flagged
                if (chars_seen[to_int(edge)] == 1) {
                    edge = tolower(edge);
                }
            } else {
                sdsl::util::set_to_value(chars_seen, 0);
            }

            W_str += edge;

            // Update loop variables
            last_label = label;
            chars_seen[to_int(edge)] = 1;
        }

        // F -> vector of integers (copied from csa)
        for (int i = 0; i < csa.sigma-1; i++) {
            F_iv[i] = ((csa.C[i+1] - k) < this->m) ? (csa.C[i+1] - k) : 0;
        }

        sdsl::util::bit_compress(F_iv);
        this->F = F_iv;

        sdsl::util::clear(csa);

        // L -> stored into an rrr_vector (plus support for rank and select)
        sdsl::util::assign(this->L, sdsl::rrr_vector<>(L_bv));
        sdsl::util::clear(L_bv);
        this->initialize_support(); 

        this->n = L_rank(this->m, 1);                                

        // W -> stored into a wavelet tree (already supporting rank and select)
        sdsl::wt_hutu<sdsl::rrr_vector<63>> W;
        sdsl::construct_im(W, W_str, 1);

        this->W = W;
    }

    // Initializes rank and select over L
    void initialize_support() {
        this->m_L_rank_0 = new sdsl::rrr_vector<>::rank_0_type (&(this->L));
        this->m_L_rank_1 = new sdsl::rrr_vector<>::rank_1_type (&(this->L));
        
        this->m_L_select_0 = new sdsl::rrr_vector<>::select_0_type (&(this->L));
        this->m_L_select_1 = new sdsl::rrr_vector<>::select_1_type (&(this->L));
    }

    // Returns the number of occurrences of c in L[0..i) 
    int L_rank(int i, int c) {
        if (i >= 0) {
            if (c == 0) return (*this->m_L_rank_0)(i);
            if (c == 1) return (*this->m_L_rank_1)(i);
        }

        return -1;
    }

    // Returns the position of the i-th occurrence of c in L
    int L_select(int i, int c) {
        if (i > 0) {
            if (c == 0) return (*this->m_L_select_0)(i);
            if (c == 1) return (*this->m_L_select_1)(i);
        }

        return -1;
    }

    // Returns the number of occurrences of c in W[0..i)
    int W_rank(int i, char c) {
        return this->W.rank(i, c);
    }

    // Returns the position of the i-th occurrence of c in L
    int W_select(int i, char c) {
        return this->W.select(i, c);
    }

    // Given an edge index i in [0, m), returns the last character of the edge and the first occurrence of that character
    std::tuple<char, int> F_lookup(int i) {
        char last_char;
        int fo = -1;

        if (this->F[0] <= i && i < this->F[1]) {
            last_char = '$';
            fo = this->F[0];
        } else if (this->F[1] <= i && i < this->F[2]) {
            last_char = 'A';
            fo = this->F[1];
        } else if (this->F[2] <= i && i < this->F[3]) {
            last_char = 'C';
            fo = this->F[2];
        } else if (this->F[3] <= i && i < this->F[4]) {
            last_char = 'G';
            fo = this->F[3];
        } else if (this->F[4] <= i) {
            last_char = 'T';
            fo = this->F[4];
        } 

        return {last_char, fo};
    }

    // Given a node index v in [0, n), returns the index i in [0, m) of the exiting edge such that L[i] == 1
    int edge_id(int v) { 
        return L_select(v+1, 1); 
    }

    // Given an edge index i in [0, m) exiting node v, returns the corresponding node index v
    int node_id(int i) { 
        return L_rank(i, 1); 
    }

    // Returns the index of the last edge of the node pointed to by edge i (returns 0 if i == p)
    int forward(int i) {
        int pointed = -1;

        if (i >= 0 && i < this->m && i != m_p) {
            char edge = toupper(this->W[i]);
            int r = W_rank(i+1, edge);
            int fo = this->F[to_int(edge)];
                    
            pointed = L_select(L_rank(fo, 1) + r, 1);

            // Special care needed to take into account the terminator $_M that does not point to anything
            if (edge == '$' && i < m_p) pointed++;
        }
        
        return pointed;
    }

    // Returns the index of the first edge that points to the node that the edge at i exits (returns p if i == 0)
    int backward(int i) {
        char last_char;
        int fo, pointed = -1;

        if (i > 0 && i < this->m) {
            std::tie (last_char, fo) = F_lookup(i);
            int j = L_rank(i, 1) - L_rank(fo, 1) + 1;        

            pointed = W_select(j, last_char);

            // Special care needed to take into account the terminator that does not point to anything
            if (last_char == '$' && pointed <= m_p) pointed = W_select(j-1, last_char);
        }

        return pointed;
    }

    // Returns the number of outgoing edges from node v
    int outdegree(int v) {
        return L_select(v+1, 1) - L_select(v, 1);
    }

    // Returns the target node after traversing edge c from node v, which might not exist
    int outgoing(int v, char c) {
        int i = edge_id(v), out = -1;

        // Check for non flagged character
        int pos_last_char_upper = W_select(W_rank(i+1, toupper(c)), toupper(c));
        // Check for flagged character
        int pos_last_char_lower = W_select(W_rank(i+1, tolower(c)), tolower(c));

        if (pos_last_char_upper > edge_id(v-1) && pos_last_char_upper <= i) {
            out = node_id(forward(pos_last_char_upper));
        } else if (pos_last_char_lower > edge_id(v-1) && pos_last_char_lower <= i) {
            out = node_id(forward(pos_last_char_lower));
        }
        
        return out;
    }

    // Returns the label of node v
    std::string label(int v) {
        std::string node_label;

        int fo, i = edge_id(v);
        char last_char;

        for (int j = 0; j < k; j++) {
            if (i == -1) {
                last_char = '$';
            } else {
                std::tie (last_char, fo) = F_lookup(i);
                i = backward(i);
            }
            
            node_label = last_char + node_label;
        }

        return node_label;
    }

    // Return the number of incoming edges
    int indegree(int v) {
        int i = edge_id(v), indeg = -1;

        int first_incoming_edge_id = backward(i);

        if (first_incoming_edge_id != -1) {
            char incoming_edge = this->W[first_incoming_edge_id];
            int num_of_edges_prev = W_rank(first_incoming_edge_id+1, incoming_edge);
        
            int next_edge;
            if (num_of_edges_prev == W_rank(this->m, incoming_edge)) {
                next_edge = m;
            } else {
                next_edge = W_select(num_of_edges_prev+1, incoming_edge);
            }

            indeg = W_rank(next_edge, tolower(incoming_edge)) - W_rank(first_incoming_edge_id, tolower(incoming_edge)) + 1;
        }

        return indeg;
    }

    // Returns the predecessor node that begins with the provided symbol
    int incoming(int v, char c) {
        int i = edge_id(v), _indegree = indegree(v), inc = -1;
        int first_incoming_edge = backward(i);

        if (first_incoming_edge != -1) {
            char edge = this->W[first_incoming_edge];   // already non flagged because of backward

            // Check first incoming edge (which has a non flagged edge)
            int candidate_node = node_id(first_incoming_edge);
            if (label(candidate_node).substr(0, 1) == std::string() + c) {
                    return candidate_node;
            }
            
            // Check following flagged edges
            int num_of_edges_prev = W_rank(first_incoming_edge, tolower(edge));

            for (int j = 1; j < _indegree; j++) { 
                candidate_node = node_id(W_select(num_of_edges_prev+j, tolower(edge)));

                if (label(candidate_node).substr(0, 1) == std::string() + c) {
                    return candidate_node;
                }
            }
        }
        
        return inc;
    }

    // Returns the index of the node labelled with string s
    int index(std::string s) {
        
        // First character of s
        char c = *(s.substr(0, 1)).c_str();

        // Range of node labels ending with s[0]
        int l = this->F[to_int(c)]; 
        int r = c != 'T' ? this->F[to_int(c) + 1] - 1 : m - 1;

        for (int i = 1; i < k; i++) {
            c = *(s.substr(i, 1)).c_str();

            // First and last occurences of c in the previous range
            int first_c = min_in_range(
                W_select(W_rank(l, c) + 1, c),
                W_select(W_rank(l, tolower(c)) + 1, tolower(c)),
                l, 
                r
            );

            int last_c = max_in_range(
                W_select(W_rank(r + 1, c), c), 
                W_select(W_rank(r + 1, tolower(c)), tolower(c)),
                l,
                r
            );

            // If there are no outgoing edges labelled with c
            if (first_c == -1 || last_c == -1) return -1;

            // New range
            l = edge_id(outgoing(node_id(first_c), c));
            r = edge_id(outgoing(node_id(last_c), c));
        }

        if (node_id(l) != node_id(r)) return -1;

        return node_id(l);
    }

    // Calculates total memory usage 
    int memory_usage_in_bits() {
        int mem_bytes = size_in_bytes(this->W) +
                        size_in_bytes(this->F) +
                        size_in_bytes(this->L);

        return 8*mem_bytes + 3*32;
    }

    // Prints all informations on the graph
    void print_graph(bool all_info=false) {
        std::cout << "k: " << this->k << "    m: " << this->m << "    n: " << this->n << '\n';

        if (all_info) {
            std::cout << "F " << this->F << std::endl << '\n';
            std::cout << "v     i     L    label     W    forward    backward   outdegree  outgoing (c=A)   indegree    incoming (c=A)" << '\n';

            for (int i = 0; i < this->m; i++) {
                std::cout.width(6);
                std::cout << std::left << node_id(i);
                
                std::cout.width(6);
                std::cout << std::left << i;

                std::cout.width(6);
                std::cout << std::left << this->L[i];

                std::cout.width(9);
                std::cout << std::left << label(node_id(i));

                std::cout.width(8);
                std::cout << std::left << this->W[i];

                std::cout.width(11);
                std::cout << std::left << forward(i);

                std::cout.width(12);
                std::cout << std::left << backward(i);

                std::cout.width(11);
                std::cout << std::left << this->outdegree(node_id(i));

                std::cout.width(16);
                std::cout << std::left << outgoing(node_id(i), '$');

                std::cout.width(10);
                std::cout << std::left << this->indegree(node_id(i));

                std::cout.width(16);
                std::cout << std::left << incoming(node_id(i), '$');

                std::cout << '\n';
            }

            std::cout  << '\n';
        }

        std::cout << "Memory used in bits: " << memory_usage_in_bits() << std::endl;
    }

};




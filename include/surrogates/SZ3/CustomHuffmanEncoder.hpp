#ifndef _HUFFMAN_ENCODER_HPP
#define _HUFFMAN_ENCODER_HPP

#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"

namespace SZ3 {

template<class T>
class CustomHuffmanEncoder : public concepts::EncoderInterface<T> {
public:

    typedef struct node_t {
        struct node_t *left, *right;
        size_t freq;
        char t; //in_node:0; otherwise:1
        T c;
    } *node;

    typedef struct HuffmanTree {
        unsigned int stateNum;
        unsigned int allNodes;
        struct node_t *pool;
        node *qqq, *qq; //the root node of the HuffmanTree is qq[1]
        int n_nodes; //n_nodes is for compression
        int qend;
        uint64_t **code;
        unsigned char *cout;
        int n_inode; //n_inode is for decompression
        int maxBitCount;
    } HuffmanTree;

    CustomHuffmanEncoder();
    ~CustomHuffmanEncoder();

    void preprocess_encode(const std::vector<T> &bins, int stateNum);
    void preprocess_encode(const T *bins, size_t num_bin, int stateNum);
    void save(unsigned char *&c);
    size_t size_est();
    size_t encode(const std::vector<T> &bins, unsigned char *&bytes);
    size_t encode(const T *bins, size_t num_bin, unsigned char *&bytes);
    void postprocess_encode();
    void preprocess_decode();
    std::vector<T> decode(const unsigned char *&bytes, size_t targetLength);
    void postprocess_decode();
    void load(const unsigned char *&c, size_t &remaining_length);
    bool isLoaded();
    unsigned char get_cout_for_state(T bin_index);

private:

    HuffmanTree *huffmanTree = NULL;
    node_t *treeRoot;
    unsigned int nodeCount;
    unsigned char sysEndianType;
    bool loaded;
    T offset;

    // Private method declarations
    HuffmanTree *createHuffmanTree(int stateNum);
    node_t *reconstruct_HuffTree_from_bytes_anyStates(const unsigned char *bytes, uint32_t nodeCount);
    node_t *new_node(size_t freq, T c, node_t *a, node_t *b);
    node_t *new_node2(T c, unsigned char t);
    void qinsert(node_t *n);
    node_t *qremove();
    void build_code(node_t *n, int len, uint64_t out1, uint64_t out2);
    void init(const T *s, size_t length);
    template<class T1> void pad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node_t *root);
    template<class T1> void unpad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node_t *root);
    template<class T1> unsigned int convert_HuffTree_to_bytes_anyStates(unsigned int nodeCount, unsigned char *out);
    void SZ_FreeHuffman();
};

} // namespace SZ3

#endif

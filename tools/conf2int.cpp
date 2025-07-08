#include <iostream>

#include "rocci_config.hpp"

ALGO algo = ALGO_NULL;
uint8_t dataType = ROCCI_FLOAT;
double eb;

signed main() {

    printf("Your Compressor: NULL\n");
    printf("Your Data Type: FLOAT\n");
    printf("Your Absolute Error Bound: ");
    scanf("%lf", &eb);
    
    Config conf;
    conf.N = 1;
    conf.dims.clear();
    // auto dims = {500, 500, 100};
    // conf.setDims(dims.begin(), dims.end());
    // conf.cmprAlgo = algo;
    conf.dataType = dataType;
    conf.absErrorBound = eb;

    uint a[10] = {0};
    unsigned char *c = reinterpret_cast<unsigned char *>(&a[0]);
    conf.save(c);

    for (int i = 0; i < 10; i++) {
        printf("%u ", a[i]);
    }
}
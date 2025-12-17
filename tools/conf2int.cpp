#include <iostream>

#include "rocci_config.hpp"

int method_id;
double eb;

signed main() {

    printf("Type Your Method Number: \n");

    printf("1  fix error bound, looking for best psnr. \n");
    printf("3* [recommend] fix error bound, looking for best compression ratio. \n");
    printf("7  fix psnr, looking for best compression ratio. \n");
    printf("12 fix compresison ratio, looking for best error bound. \n");
    printf("13 fix compression ratio, looking for best psnr. \n");

    scanf("%d", &method_id);

    printf("Your Compression Parameter (either error bound, psnr or compresison ratio): \n");
    scanf("%lf", &eb);
    
    Config conf;
    conf.N = 1;
    conf.dims.clear();
    // auto dims = {500, 500, 100};
    // conf.setDims(dims.begin(), dims.end());
    // conf.cmprAlgo = algo;
    // conf.dataType = dataType;
    conf.method_id = method_id;
    conf.absErrorBound = eb;

    uint a[10] = {0};
    unsigned char *c = reinterpret_cast<unsigned char *>(&a[0]);
    conf.save(c);

    printf("-UD=32088,0,3");
    for (int i = 0; i < 3; i++) {
        printf(",%u", a[i]);
    }
    printf("\n");

    return 0;
}
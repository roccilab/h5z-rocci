//
// Created by Kai Zhao on 12/9/19.
//

#include "rocci_FileUtil.h"

void readfile(const char *file, const size_t num, size_t elem_size, void *data) {
    FILE *fin = fopen(file, "rb");
    if (!fin) {
        printf("Error, Couldn't find the file: %s\n", file);
        exit(0);
    }
    fseek(fin, 0, SEEK_END);
    const size_t num_elements = ftell(fin) / elem_size;
    assert(num_elements == num && "File size is not equal to the input setting");
    fseek(fin, 0, SEEK_SET);
    fread(data, elem_size, num_elements, fin);
    fclose(fin);
}


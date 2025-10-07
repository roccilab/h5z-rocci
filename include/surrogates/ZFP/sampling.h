#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <type_traits>
#include <cmath>
#include <numeric>
#include <stdio.h>
#include <chrono>
#include <cmath>
#include <cfloat>
#include <sys/stat.h>

/*
This is a block sampling algorithm, I used a derivation of an algorithm 
used to interpolate colors in graphics to allow us to get more blocks, i.e.
get closer to the target sample rate. But in principle uniform sampling
is fine too. For ZFP we just want blocks of size 4^n so this conceptually attempts to
sample blocks rather than individual data points which is why the indices in the loops
are the way that they are - the inner most loop is just copying whatever blocks we ended up sampling
to the output buffer for the sampled dataset.
*/
// 4^n block sampling
template<class T>
std::vector<T>
sampling_1d(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, float sample_ratio) {
	size_t nbEle = dims[0];

	// get num elements to sample sample_ratio, floor to nearest multiple of 4
	size_t temp = nbEle * sample_ratio;
	sample_num = temp - (temp % 4);
	sample_dims.push_back(sample_num);

	size_t nbBlocks = sample_num / 4;

	// compute stride for blocks
	size_t totalBlocks = nbEle/4;
	// totalBlocks = (totalBlocks - (totalBlocks%4)) / 4;
	size_t stride_0 = roundf(((float)totalBlocks / (float)nbBlocks)) * 4;

	std::vector<T> sample_data(sample_num, 0);
	size_t idx = 0;
	for(size_t i = 0; i + 4 < nbEle;){

		for(int j = 0; j < 4; j++)
			sample_data[idx++] = data[i+j];

		i += stride_0;

		if(idx == sample_num) break;

	}

	return sample_data;

}

template<class T>
std::vector<T>
sampling_2d(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, float sample_ratio) {
	size_t nbEle = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
	// compute number of blocks to sample in order to roughly maintain dimension ratio
	size_t blocksX = dims[0] / 4;
	size_t blocksY = dims[1] / 4;
	size_t nbBlocksX = sqrt(sample_ratio)*blocksX;
	size_t nbBlocksY = sqrt(sample_ratio)*blocksY;

	sample_num = nbBlocksX*nbBlocksY*16;
	sample_dims.push_back(nbBlocksX*4);
	sample_dims.push_back(nbBlocksY*4);

	std::vector<T> sample_data(sample_num, 0);

	for(size_t i = 0; i < nbBlocksX; i++){
		size_t x_block_ind = roundf(((float)i/(nbBlocksX-1)) * (blocksX-1));
		size_t x_ind = x_block_ind * 4;
		for(size_t j = 0; j < nbBlocksY; j++){
			size_t y_block_ind = roundf(((float)j/(nbBlocksY-1)) * (blocksY-1));
			size_t y_ind = y_block_ind * 4;

			for(size_t ii = 0; ii < 4; ii++) {
				for(size_t jj = 0; jj < 4; jj++) {
					size_t data_idx = (x_ind+ii)*dims[1] + y_ind + jj;
					size_t scal_i = i * 4;
					size_t scal_j = j * 4;
					size_t sample_idx = (scal_i+ii)*sample_dims[1] + scal_j + jj;

					sample_data[sample_idx] = data[data_idx];
				}
			}
		}
	}

	return sample_data;

}	


template<class T>
std::vector<T>
sampling_3d(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, float sample_ratio) {
	size_t nbEle = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
	// compute number of blocks to sample in order to roughly maintain dimension ratio
	size_t blocksX = dims[0] / 4;
	size_t blocksY = dims[1] / 4;
	size_t blocksZ = dims[2] / 4;
	size_t nbBlocksX = std::cbrt(sample_ratio)*blocksX;
	size_t nbBlocksY = std::cbrt(sample_ratio)*blocksY;
	size_t nbBlocksZ = std::cbrt(sample_ratio)*blocksZ;

	sample_num = nbBlocksX*nbBlocksY*nbBlocksZ*64;
	sample_dims.push_back(nbBlocksX*4);
	sample_dims.push_back(nbBlocksY*4);
	sample_dims.push_back(nbBlocksZ*4);


	std::vector<T> sample_data(sample_num, 0);

	for(size_t i = 0; i < nbBlocksX; i++){
		size_t x_block_ind = roundf(((float)i/(nbBlocksX-1)) * (blocksX-1));
		size_t x_ind = x_block_ind * 4;

		for(size_t j = 0; j < nbBlocksY; j++){
			size_t y_block_ind = roundf(((float)j/(nbBlocksY-1)) * (blocksY-1));
			size_t y_ind = y_block_ind * 4;

			for(size_t k = 0; k < nbBlocksZ; k++){
				size_t z_block_ind = roundf(((float)k/(nbBlocksZ-1)) * (blocksZ-1));
				size_t z_ind = z_block_ind * 4;

				for(size_t ii = 0; ii < 4; ii++) {
					for(size_t jj = 0; jj < 4; jj++) {
						for(size_t kk = 0; kk < 4; kk++){
							size_t data_idx = (x_ind+ii)*dims[1]*dims[2] + (y_ind+jj)*dims[2] + z_ind + kk;
							size_t scal_i = i * 4;
							size_t scal_j = j * 4;
							size_t scal_k = k * 4;
							size_t sample_idx = 
								(scal_i+ii)*sample_dims[1]*sample_dims[2] + (scal_j+jj)*sample_dims[2] + scal_k + kk;

							sample_data[sample_idx] = data[data_idx];
						}
					}
				}

			}
			
		}
	}

	return sample_data;

}	

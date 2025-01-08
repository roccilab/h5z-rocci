#include <ZFP/estimate_metric.h>
#include <ZFP/sampling.h>
#include <qcat_ssim.h>
#include <iostream>
#include <math.h>


double zfp_estimate_cr_float(const struct pressio_data* input, double eb, double sample_ratio) {
    assert(input->num_dimensions() < 4);
    std::vector<size_t> dims = input->dimensions();
      size_t sample_num;
      std::vector<size_t> sample_dims;
      std::vector<float> sample_data;
      float * data_ = static_cast<float*>(input->data());
      switch(input->num_dimensions()){
			case 1:
				sample_data = sampling_1d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			case 2:
				sample_data = sampling_2d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			case 3:
				sample_data = sampling_3d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			default:
				std::cerr << "Dims must be between size 1 and 3" << std::endl;
                exit(1);
		}
      
      float* sampled_dataset = sample_data.data();
      struct pressio_data* input_data =
				pressio_data_new_move(pressio_float_dtype, sampled_dataset, sample_dims.size(), sample_dims.data(),
          pressio_data_libc_free_fn, NULL);
      struct pressio_data* compressed_data =
    			pressio_data_new_empty(pressio_byte_dtype, 0, NULL);

        pressio_compressor comp = compressor_plugins().build("zfp");
      struct pressio_options* cmp_options = pressio_compressor_get_options(&comp);
      pressio_options_set_double(cmp_options, "zfp:accuracy", eb);
      pressio_compressor_set_options(&comp, cmp_options);
      if(comp->compress(input_data, compressed_data)) {
            std::cerr << comp->error_msg() << std::endl;
            exit(1);
        }
      double estCR = (double)input_data->size_in_bytes()/compressed_data->size_in_bytes();
      return estCR;
    
}

double calc_psnr(float* ori, float* other, size_t nbEle){
    float mse = 0.0;
    float eps = 1e-16;
    float maxVal = std::numeric_limits<float>::min();
    float minVal = std::numeric_limits<float>::max();
    float mean_err = 0.0;
    for(int i = 0; i < nbEle; i++){
        mse += pow(ori[i] - other[i], 2);
        mean_err += abs(ori[i] - other[i]);
        if(ori[i] > maxVal) maxVal = ori[i];
        if(ori[i] < minVal) minVal = ori[i];
    }
    // printf("sqerr=%f\n", mse);
    mse /= nbEle;
    mean_err /= nbEle;
    float value_range = maxVal - minVal;
    double psnr = -20.0*log10((sqrt(mse)/value_range) + eps);
    // printf("max=%f, mse=%f, psnr=%f\n", value_range, mse, psnr);
    // printf("nbEle=%i\n", nbEle);
    // printf("mean_err=%.15f\n", mean_err);
    return psnr;
}

double zfp_estimate_psnr_float(const struct pressio_data* input, double eb, double sample_ratio) {
    assert(input->num_dimensions() < 4);
    std::vector<size_t> dims = input->dimensions();
      size_t sample_num;
      std::vector<size_t> sample_dims;
      std::vector<float> sample_data;
      float * data_ = static_cast<float*>(input->data());
      switch(input->num_dimensions()){
			case 1:
				sample_data = sampling_1d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			case 2:
				sample_data = sampling_2d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			case 3:
				sample_data = sampling_3d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			default:
				std::cerr << "Dims must be between size 1 and 3" << std::endl;
                exit(1);
		}
      
      float* sampled_dataset = sample_data.data();
      struct pressio_data* input_data =
				pressio_data_new_move(pressio_float_dtype, sampled_dataset, sample_dims.size(), sample_dims.data(),
          pressio_data_libc_free_fn, NULL);
      struct pressio_data* compressed_data =
    			pressio_data_new_empty(pressio_byte_dtype, 0, NULL);

        pressio_compressor comp = compressor_plugins().build("zfp");
      struct pressio_options* cmp_options = pressio_compressor_get_options(&comp);
      pressio_options_set_double(cmp_options, "zfp:accuracy", eb);
      pressio_compressor_set_options(&comp, cmp_options);
        if(comp->compress(input_data, compressed_data)) {
            std::cerr << comp->error_msg() << std::endl;
            exit(1);
        }

        pressio_data output = pressio_data::owning(input->dtype(), input->dimensions());
        if(comp->decompress(compressed_data, &output)) {
            std::cerr << comp->error_msg() << std::endl;
            exit(1);
        }

    float* decData = static_cast<float*>(output.data());

    return calc_psnr(sampled_dataset, decData, sample_num);
}

double zfp_estimate_ssim_float(const struct pressio_data* input, double eb, double sample_ratio) {
    assert(input->num_dimensions() < 4);
    std::vector<size_t> dims = input->dimensions();
      size_t sample_num;
      std::vector<size_t> sample_dims;
      std::vector<float> sample_data;
      float * data_ = static_cast<float*>(input->data());
      switch(input->num_dimensions()){
			case 1:
				sample_data = sampling_1d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			case 2:
				sample_data = sampling_2d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			case 3:
				sample_data = sampling_3d<float>(data_, dims, sample_num, sample_dims, sample_ratio);
				break;
			default:
				std::cerr << "Dims must be between size 1 and 3" << std::endl;
                exit(1);
		}
      
      float* sampled_dataset = sample_data.data();
      struct pressio_data* input_data =
				pressio_data_new_move(pressio_float_dtype, sampled_dataset, sample_dims.size(), sample_dims.data(),
          pressio_data_libc_free_fn, NULL);
      struct pressio_data* compressed_data =
    			pressio_data_new_empty(pressio_byte_dtype, 0, NULL);


    pressio_compressor comp = compressor_plugins().build("zfp");
      struct pressio_options* cmp_options = pressio_compressor_get_options(&comp);
      pressio_options_set_double(cmp_options, "zfp:accuracy", eb);
      pressio_compressor_set_options(&comp, cmp_options);
        if(comp->compress(input_data, compressed_data)) {
            std::cerr << comp->error_msg() << std::endl;
            exit(1);
        }

        pressio_data output = pressio_data::owning(input->dtype(), input->dimensions());
        if(comp->decompress(compressed_data, &output)) {
            std::cerr << comp->error_msg() << std::endl;
            exit(1);
        }

    float* decData = static_cast<float*>(output.data());

    auto n_dims = input->num_dimensions();
    double estimatedSSIM;
    if(n_dims == 1){
        estimatedSSIM = SSIM_1d_windowed_float(sampled_dataset, decData, sample_dims[0], 8, 8);
    }
    else if (n_dims == 2){
        estimatedSSIM = SSIM_2d_windowed_float(sampled_dataset, decData, sample_dims[1], sample_dims[0], 8,8, 8,8);
    }
    else { // n_dims == 3
        estimatedSSIM = SSIM_3d_windowed_float(sampled_dataset, decData, sample_dims[2], sample_dims[1], sample_dims[0], 8,8,8, 8,8,8);
    }

    return estimatedSSIM;
}

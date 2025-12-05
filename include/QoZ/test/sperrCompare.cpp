#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <cfloat>
#include <cmath>
#include "QoZ/api/sz.hpp"

using namespace QoZ;

// params for experiments
size_t blocksizes[] = {16,32,64};
bool skip_setting[] = {false, true};
double sample_rates[] = {0.01, 0.10, 0.20, 0.30, 0.40};

template<class T, uint N>
std::vector<double> estimate_compress(Config conf, T *data, double abs, std::string file, std::ofstream& output, std::ofstream& sampleout){
	conf.errorBoundMode = QoZ::EB_ABS;
    conf.absErrorBound = abs;

	std::vector<double> results;

	char ebstr[2048];
	std::sprintf(ebstr, "%0.10f", abs);

	size_t outSize = 0;
	printf("startcmpr\n");
	Timer timer(true);
	auto cmpData = SPERR_Compress<T, N>(conf, data, outSize);
	double dur = timer.stop();

	printf("cmpr time: %f s\n\n", dur);

	double trueCR = conf.num * sizeof(T) * 1.0 / outSize;

	output << file << "," << ebstr << "," << std::to_string(trueCR) << "," << std::to_string(dur) << "\n" << std::flush;

	delete[] cmpData;

	// results.push_back(trueCR);
	// results.push_back(dur);

	for(double sample_rate : sample_rates){
		printf("SR: %f\n", sample_rate);
		for(bool skip_outlier : skip_setting){
			printf("Skip: %i\n", skip_outlier);
			for(size_t bs : blocksizes){
				printf("BS: %i\n", bs);
				timer.start();
				T* data_cpy = new T[conf.num];
				memcpy(data_cpy, data, sizeof(T)*conf.num);
				double cpydur = timer.stop("time to copy");

				printf("copy time: %f s\n", cpydur);

				timer.start();
				double estCR = estimateSPERRCRbasedonErrorBound<T,N>(abs, data_cpy, sample_rate, bs, conf.dims, skip_outlier);
				double estDur = timer.stop();

				printf("est time: %f s\n\n", estDur);

				// timer.start();
				// results.push_back(estCR);
				// timer.stop("push to estCR");
				// timer.start();
				// results.push_back(estDur);
				// timer.stop("push to estDUR");

				char ratestr[256];
				sprintf(ratestr, "%f", sample_rate);
				char skipstr[256];
				std::sprintf(skipstr, "%i", skip_outlier);
				char blockstr[256];
				std::sprintf(blockstr, "%i", bs);

				sampleout << file << "," << ebstr << "," << ratestr << "," << skipstr << "," << blockstr << "," << std::to_string(estCR) << "," << std::to_string(estDur) << "\n" << std::flush;




				delete[] data_cpy;
			}
		}
	}

	return results;

}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cout << "usage: " << argv[0] << " data_file r1 r2 r3 tag" << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat 33120 69 69 qmcpack" << std::endl;
        return 0;
    }

	std::string targetFiles = argv[1];
	std::ifstream file(targetFiles);
    std::string str;
    std::vector<std::string> files;
    while (std::getline(file, str))
    {
        files.push_back(str);
    }
    file.close();

	std::string tag;
	size_t r1,r2,r3,r4,r5;
	r1 = atoi(argv[2]);
	r2 = atoi(argv[3]);
	r3 = atoi(argv[4]);
	tag = argv[5];

	size_t dim;
	if(r1 == 0){
		printf("Can't have 0-dim data\n");
		return 1;
	}
	else if(r2 == 0){
		dim = 1;
	}
	else if(r3 == 0){
		dim = 2;
	} 
	else {
		dim = 3;
	}


	printf("Dim: %i", dim);
	// create dir if not exists
	std::string path = "./analysis";
	struct stat sb;
	if(stat(path.c_str(), &sb) != 0){
		mkdir(path.c_str(), 0777);
	}
	std::string outpath = "./analysis/" + tag + "/";

	struct stat sb2;
	if(stat(outpath.c_str(), &sb2) != 0){
		mkdir(outpath.c_str(), 0777);
	}

	// std::string outfile = outpath + tag + ".csv";
	// std::ofstream output(outfile);

	std::string outfile = outpath + tag + ".csv";
	std::ofstream output(outfile);
	
	std::string samplefile = outpath + tag + "_estimate.csv";
	std::ofstream sampleout(samplefile);

	std::vector<float> eb;
    for(int i = 2; i < 6; i++){
        eb.push_back(pow(10, i+1)*10e-9);
        // eb.push_back(pow(10, i+1)*1);
    }

    std::cout << eb[0] << " " << eb[1] << " " << eb[2] << " ... " << eb[eb.size()-1] << std::endl;
	for(int i = 0; i < files.size(); i++){
		if(i >= 100 ) break;
        std::cout << files[i] << std::endl;
        size_t num = 0;
        auto data = QoZ::readfile<float>(files[i].c_str(), num);
        float min = FLT_MAX;
        float max = FLT_MIN;
        for(int i = 0; i < num; i++){
            if(data[i] < min) min = data[i];
            if(data[i] > max) max = data[i];
        }
        float range = fabs(max - min);
        printf("Data Range: %f\n", range);

		for(int j = 0; j < eb.size(); j++){
			printf("inner loop\n");
			float* data_ = new float[sizeof(float)*num];
            memcpy(data_, data.get(), sizeof(float)*num);

			float eb_ = eb[j] * range;
            printf("EB: %.10f\n", eb_);

			std::vector<double> results;
			printf("branch\n");
			if (dim == 1) {
				Config config(r1);
				results = estimate_compress<float, 1>(config, data_, eb_, files[i], output, sampleout);
			} else if (dim == 2) {
				Config config(r2, r1);
				results = estimate_compress<float, 2>(config, data_, eb_, files[i], output, sampleout);
			} else if (dim == 3) {
				Config config(r3, r2, r1);
				results = estimate_compress<float, 3>(config, data_, eb_, files[i], output, sampleout);
			}

			printf("exit branch\n");


			// printf("write\n");
            // char ebstr[2048];
            // std::sprintf(ebstr, "%0.10f", eb_);//eb[j]);

			// output << files[i] << "," << ebstr << "," << std::to_string(results[0]) << "," << std::to_string(results[1]) << "\n" << std::flush;
			// int k = 2;
			// for(double rate : sample_rates){
			// 	char ratestr[256];
			// 	sprintf(ratestr, "%f", rate);
			// 	for(bool skip : skip_setting){
			// 		char skipstr[256];
			// 		std::sprintf(skipstr, "%i", skip);
			// 		for(size_t bs : blocksizes){
			// 			char blockstr[256];
			// 			std::sprintf(blockstr, "%i", bs);

			// 			sampleout << files[i] << "," << ebstr << "," << ratestr << "," << skipstr << "," << blockstr << "," << std::to_string(results[k]) << "," << std::to_string(results[k+1]) << "\n" << std::flush;
			// 			k += 2;
			// 		}
			// 	}
			// }
			// printf("exit write\n");

            delete[] data_;

		}

	}

    // output.close();
    output.close();
    sampleout.close();



}
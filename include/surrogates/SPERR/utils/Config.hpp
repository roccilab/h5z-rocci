//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_Config_HPP
#define SZ_Config_HPP

#include <iostream>
#include <vector>
#include <numeric>
#include "QoZ/def.hpp"
#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/utils/inih/INIReader.h"
#include <SPERR/sperr/Conditioner.h>

namespace QoZ {

    enum EB {
        EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL
    };
    constexpr const char *EB_STR[] = {"ABS", "REL", "PSNR", "NORM", "ABS_AND_REL", "ABS_OR_REL"};
    constexpr EB EB_OPTIONS[] = {EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL};

    enum ALGO {
        ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP,ALGO_INTERP_BLOCKED
    };
    const char *ALGO_STR[] = {"ALGO_LORENZO_REG", "ALGO_INTERP_LORENZO", "ALGO_INTERP","ALGO_INTERP_BLOCKED"};
    constexpr const ALGO ALGO_OPTIONS[] = {ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_INTERP_BLOCKED};

    enum INTERP_ALGO {
        INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC
    };
    constexpr const char *INTERP_ALGO_STR[] = {"INTERP_ALGO_LINEAR", "INTERP_ALGO_CUBIC"};
    constexpr const INTERP_ALGO INTERP_ALGO_OPTIONS[] = { INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC };

    enum TUNING_TARGET {
        TUNING_TARGET_RD, TUNING_TARGET_CR, TUNING_TARGET_SSIM, TUNING_TARGET_AC
    };
    constexpr const char *TUNING_TARGET_STR[] = {"TUNING_TARGET_RD", "TUNING_TARGET_CR", "TUNING_TARGET_SSIM", "TUNING_TARGET_AC"};
    constexpr const TUNING_TARGET TUNING_TARGET_OPTIONS[] = {
        TUNING_TARGET_RD, TUNING_TARGET_CR, TUNING_TARGET_SSIM, TUNING_TARGET_AC
    };

    template<class T>
    const char *enum2Str(T e) {
        if (std::is_same<T, ALGO>::value) {
            return ALGO_STR[e];
        } else if (std::is_same<T, INTERP_ALGO>::value) {
            return INTERP_ALGO_STR[e];
        } else if (std::is_same<T, EB>::value) {
            return EB_STR[e];

        }  else if (std::is_same<T, TUNING_TARGET>::value) {
            return TUNING_TARGET_STR[e];
            
        }else {
            printf("invalid enum type for enum2Str()\n ");
            exit(0);
        }
    }

    class Config {
    public:
        template<class ... Dims>
        Config(Dims ... args) {
            dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
            coeffs_dims=dims;
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
            coeffs_num = num;
            blockSize = (N == 1 ? 128 : (N == 2 ? 16 : 6));
            pred_dim = N;
            stride = blockSize;
        }

        template<class Iter>
        size_t setDims(Iter begin, Iter end) {
            dims = std::vector<size_t>(begin, end);
            coeffs_dims=dims;
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
            coeffs_num = num;
            //added
            blockSize = (N == 1 ? 128 : (N == 2 ? 16 : 6));
            pred_dim = N;
            stride = blockSize;
            //added end 
            return num;
        }

        void loadcfg(const std::string &cfgpath) {
            INIReader cfg(cfgpath);

            if (cfg.ParseError() != 0) {
                std::cout << "Can't load cfg file  <<" << cfgpath << std::endl;
                exit(0);
            } else {
                //std::cout << "Load cfg from " << cfgpath << std::endl;
            }

            auto cmprAlgoStr = cfg.Get("GlobalSettings", "CmprAlgo", "");
            if (cmprAlgoStr == ALGO_STR[ALGO_LORENZO_REG]) {
                cmprAlgo = ALGO_LORENZO_REG;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP_LORENZO]) {
                cmprAlgo = ALGO_INTERP_LORENZO;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP]) {
                cmprAlgo = ALGO_INTERP;
            }
            else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP_BLOCKED]) {
                cmprAlgo = ALGO_INTERP_BLOCKED;
            }
            auto ebModeStr = cfg.Get("GlobalSettings", "ErrorBoundMode", "");
            if (ebModeStr == EB_STR[EB_ABS]) {
                errorBoundMode = EB_ABS;
            } else if (ebModeStr == EB_STR[EB_REL]) {
                errorBoundMode = EB_REL;
            } else if (ebModeStr == EB_STR[EB_PSNR]) {
                errorBoundMode = EB_PSNR;
            } else if (ebModeStr == EB_STR[EB_L2NORM]) {
                errorBoundMode = EB_L2NORM;
            } else if (ebModeStr == EB_STR[EB_ABS_AND_REL]) {
                errorBoundMode = EB_ABS_AND_REL;
            } else if (ebModeStr == EB_STR[EB_ABS_OR_REL]) {
                errorBoundMode = EB_ABS_OR_REL;
            }
            auto tuningTargetStr = cfg.Get("GlobalSettings", "tuningTarget", "");
            if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_RD]) {
                tuningTarget = TUNING_TARGET_RD;
            } else if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_CR]) {
                tuningTarget = TUNING_TARGET_CR;
            }
            else if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_SSIM]) {
                tuningTarget = TUNING_TARGET_SSIM;
            }
            else if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_AC]) {
                tuningTarget = TUNING_TARGET_AC;
            }
                


            absErrorBound = cfg.GetReal("GlobalSettings", "AbsErrorBound", absErrorBound);
            relErrorBound = cfg.GetReal("GlobalSettings", "RelErrorBound", relErrorBound);
            psnrErrorBound = cfg.GetReal("GlobalSettings", "PSNRErrorBound", psnrErrorBound);
            l2normErrorBound = cfg.GetReal("GlobalSettings", "L2NormErrorBound", l2normErrorBound);
            //prewave_absErrorBound= cfg.GetReal("GlobalSettings", "prewave_absErrorBound", prewave_absErrorBound);
            alpha = cfg.GetReal("AlgoSettings", "alpha", alpha);
            beta = cfg.GetReal("AlgoSettings", "beta", beta);
            autoTuningRate = cfg.GetReal("AlgoSettings", "autoTuningRate", autoTuningRate);
            predictorTuningRate = cfg.GetReal("AlgoSettings", "predictorTuningRate", predictorTuningRate);
            wavelet_rel_coeff = cfg.GetReal("AlgoSettings", "wavelet_rel_coeff", wavelet_rel_coeff);
            preTrim = cfg.GetReal("AlgoSettings", "preTrim", preTrim);
            waveletTuningRate = cfg.GetReal("AlgoSettings", "waveletTuningRate", waveletTuningRate);
            waveletBrFix = cfg.GetReal("AlgoSettings", "waveletBrFix", waveletBrFix);
            waveletBrFix2 = cfg.GetReal("AlgoSettings", "waveletBrFix2", waveletBrFix2);
            waveletMseFix = cfg.GetReal("AlgoSettings", "waveletMseFix", waveletMseFix);
            waveletMseFix2 = cfg.GetReal("AlgoSettings", "waveletMseFix2", waveletMseFix2);
            lorenzoBrFix = cfg.GetReal("AlgoSettings", "lorenzoBrFix", lorenzoBrFix);
            //anchorThreshold= cfg.GetReal("AlgoSettings", "anchorThreshold", anchorThreshold);
            //sperr_eb_coeff = cfg.GetReal("AlgoSettings", "sperr_eb_coeff", sperr_eb_coeff);

            openmp = cfg.GetBoolean("GlobalSettings", "OpenMP", openmp);
            lorenzo = cfg.GetBoolean("AlgoSettings", "Lorenzo", lorenzo);
            lorenzo2 = cfg.GetBoolean("AlgoSettings", "Lorenzo2ndOrder", lorenzo2);
            regression = cfg.GetBoolean("AlgoSettings", "Regression", regression);
            regression2 = cfg.GetBoolean("AlgoSettings", "Regression2ndOrder", regression2);
            //external_wave = cfg.GetBoolean("AlgoSettings", "external_wave", external_wave);
            
            
            
            
            auto interpAlgoStr = cfg.Get("AlgoSettings", "InterpolationAlgo", "");
            if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_LINEAR]) {
                interpAlgo = INTERP_ALGO_LINEAR;
            } else if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_CUBIC]) {
                interpAlgo = INTERP_ALGO_CUBIC;
            }
            interpDirection = cfg.GetInteger("AlgoSettings", "InterpDirection", interpDirection);
            interpBlockSize = cfg.GetInteger("AlgoSettings", "InterpBlockSize", interpBlockSize);
            blockSize = cfg.GetInteger("AlgoSettings", "BlockSize", blockSize);
            quantbinCnt = cfg.GetInteger("AlgoSettings", "QuantizationBinTotal", quantbinCnt);
            maxStep=cfg.GetInteger("AlgoSettings", "maxStep", maxStep);
            sampleBlockSize=cfg.GetInteger("AlgoSettings", "sampleBlockSize", sampleBlockSize);
            levelwisePredictionSelection=cfg.GetInteger("AlgoSettings", "levelwisePredictionSelection", levelwisePredictionSelection);
            //exhaustiveTuning=cfg.GetInteger("AlgoSettings", "exhaustiveTuning", exhaustiveTuning);
            testLorenzo=cfg.GetInteger("AlgoSettings", "testLorenzo", testLorenzo);
            //linearReduce=cfg.GetBoolean("AlgoSettings", "linearReduce", linearReduce);
            //train=cfg.GetBoolean("AlgoSettings", "train", train);
            //useCoeff=cfg.GetBoolean("AlgoSettings", "useCoeff", useCoeff);
            //regSampleStep=cfg.GetInteger("AlgoSettings", "regSampleStep", regSampleStep);
            //multiDimInterp=cfg.GetInteger("AlgoSettings", "multiDimInterp", multiDimInterp);
            profiling=cfg.GetInteger("AlgoSettings", "profiling", profiling);
            SSIMBlockSize=cfg.GetInteger("AlgoSettings", "SSIMBlockSize", SSIMBlockSize);
            fixBlockSize=cfg.GetInteger("AlgoSettings", "fixBlockSize", fixBlockSize);
            verbose=cfg.GetBoolean("AlgoSettings", "verbose", verbose);
            QoZ=cfg.GetBoolean("AlgoSettings", "QoZ", QoZ);
            QoZ=cfg.GetBoolean("AlgoSettings", "FZ", FZ);
            sperrWithoutWave=cfg.GetBoolean("AlgoSettings", "sperrWithoutWave",sperrWithoutWave);
            pyBind=cfg.GetBoolean("AlgoSettings", "pyBind",pyBind);
            //profilingFix=cfg.GetBoolean("AlgoSettings", "profilingFix",profilingFix);



            blockwiseSampleBlockSize=cfg.GetInteger("AlgoSettings", "blockwiseSampleBlockSize", blockwiseSampleBlockSize);
            //pdTuningRealComp=cfg.GetBoolean("AlgoSettings", "pdTuningRealComp", pdTuningRealComp);
           // pdTuningAbConf=cfg.GetInteger("AlgoSettings", "pdTuningAbConf", pdTuningAbConf);
           // pdAlpha=cfg.GetReal("AlgoSettings", "pdAlpha", pdAlpha);
            //pdBeta=cfg.GetReal("AlgoSettings", "pdBeta", pdBeta);
            //lastPdTuning=cfg.GetBoolean("AlgoSettings", "lastPdTuning", lastPdTuning);
            abList=cfg.GetInteger("AlgoSettings", "abList", abList);
            crossBlock=cfg.GetInteger("AlgoSettings", "crossBlock", crossBlock);
            sampleBlockSampleBlockSize=cfg.GetInteger("AlgoSettings", "sampleBlockSampleBlockSize", sampleBlockSampleBlockSize);
            peTracking=cfg.GetBoolean("AlgoSettings", "peTracking", peTracking);
            wavelet=cfg.GetInteger("AlgoSettings", "wavelet", wavelet);
            offsetPredictor=cfg.GetInteger("AlgoSettings", "offsetPredictor", offsetPredictor);
            //transformation=cfg.GetInteger("AlgoSettings", "transformation", transformation);
            trimToZero = cfg.GetInteger("AlgoSettings", "trimToZero", trimToZero);
            pid = cfg.GetInteger("AlgoSettings", "pid", pid);
            blockOrder = cfg.GetInteger("AlgoSettings", "blockOrder", blockOrder);
            coeffTracking = cfg.GetInteger("AlgoSettings", "coeffTracking", coeffTracking);
            waveletTest = cfg.GetInteger("AlgoSettings", "waveletTest", waveletTest);
            waveletAutoTuning = cfg.GetInteger("AlgoSettings", "waveletAutoTuning", waveletAutoTuning);
            var_first = cfg.GetInteger("AlgoSettings", "var_first", var_first);
            profStride = cfg.GetInteger("AlgoSettings", "profStride", profStride);
            sperr = cfg.GetInteger("AlgoSettings", "sperr", sperr);
            waveAutoFix = cfg.GetInteger("AlgoSettings", "waveAutoFix", waveAutoFix);
            conditioning = cfg.GetInteger("AlgoSettings", "conditioning", conditioning);
            fixWave = cfg.GetInteger("AlgoSettings", "fixWave", fixWave);
           // minAnchorLevel = cfg.GetInteger("AlgoSettings", "minAnchorLevel", minAnchorLevel);




        }

        static size_t size_est() {
            
            return sizeof(size_t) * 10 + sizeof(double) * 8 + sizeof(bool) * 10 + sizeof(uint8_t) * 12 + sizeof(int) * 10 + 200; //doubled SZ3 est+100
        }

        void save(unsigned char *&c) {
            write(N, c);
            write(dims.data(), dims.size(), c);
            write(num, c);
            write(cmprAlgo, c);
            write(errorBoundMode, c);
            write(tuningTarget, c);
            write(absErrorBound, c);
            write(relErrorBound, c);
            write(alpha,c);
            write(beta,c);
            //write(autoTuningRate,c);
            //write(predictorTuningRate,c);
            write(lorenzo, c);
            write(lorenzo2, c);
            write(regression, c);
            write(regression2, c);
            write(interpAlgo, c);
            write(interpDirection, c);
            write(interpBlockSize, c);
            write(lossless, c);
            write(encoder, c);
            write(quantbinCnt, c);
            write(blockSize, c);
            
            write(levelwisePredictionSelection, c);
            write(blockwiseTuning, c);
            write(stride, c);
            write(maxStep,c);
            write(pred_dim, c);
            write(openmp, c);
            write(fixBlockSize, c);
            write(blockwiseSampleBlockSize, c);
            //write(QoZ, c);//recently changed.
            write(crossBlock, c);
            write(wavelet, c);
            write(firstSize, c);
            write(offsetPredictor, c);
            //write(transformation, c);
            //write(external_wave, c);
            
            write(coeffs_num, c);
            if(coeffs_num>0){
                size_t tempsize=coeffs_dims.size();
                write(tempsize,c);
                write(coeffs_dims.data(), tempsize, c);
            }
            write(pid, c);
            write(blockOrder, c);
            write(sperr, c);
            write(sperrWithoutWave, c);
            //write(trimToZero, c);
            //write(prewave_absErrorBound, c);
            
            write(pyBind,c);
            if(pyBind>0){
                write(metadata.size(),c);
                write(metadata.data(),metadata.size(),c);
            }
            write(conditioning, c);
            if(conditioning>0){
                meta_size=meta.size();
                //std::vector<uint8_t> temp_meta(meta.begin(),meta.end());
                write(meta_size,c);
                write(meta.data(),meta_size,c);
               
            }

            
        };

        void load(const unsigned char *&c) {
            read(N, c);
            dims.resize(N);
            read(dims.data(), N, c);
            read(num, c);
            read(cmprAlgo, c);
            read(errorBoundMode, c);
            read(tuningTarget, c);
            read(absErrorBound, c);
            read(relErrorBound, c);
            read(alpha,c);
            read(beta,c);
            //read(autoTuningRate,c);
            //read(predictorTuningRate,c);
            read(lorenzo, c);
            read(lorenzo2, c);
            read(regression, c);
            read(regression2, c);
            read(interpAlgo, c);
            read(interpDirection, c);
            read(interpBlockSize, c);
            read(lossless, c);
            read(encoder, c);
            read(quantbinCnt, c);
            read(blockSize, c);
            read(levelwisePredictionSelection, c);
            read(blockwiseTuning, c);
            read(stride, c);
            read(maxStep,c);
            read(pred_dim, c);
            read(openmp, c);
            read(fixBlockSize, c);
            read(blockwiseSampleBlockSize, c);
            //read(QoZ, c);//recently changed.
            read(crossBlock, c);
            read(wavelet, c);
            read(firstSize, c);
            read(offsetPredictor, c);
            //read(transformation, c);
            //read(external_wave, c);
            
            read(coeffs_num, c);
            if(coeffs_num>0){
                size_t tempsize;
                read(tempsize,c);
                coeffs_dims.resize(tempsize);
                read(coeffs_dims.data(), tempsize, c);
            }
            read(pid, c);
            read(blockOrder, c);
            read(sperr, c);
            read(sperrWithoutWave, c);
            //read(trimToZero, c);
            //read(prewave_absErrorBound, c);
            read(pyBind,c);
            if(pyBind>0){
             
                size_t msize;
                read(msize,c);
              
                metadata.resize(msize);
                read(metadata.data(),msize,c);
               
            }

            read(conditioning, c);
            if(conditioning>0){
               // std::cout<<"dwad"<<std::endl;
                read(meta_size,c);
                //std::cout<<meta_size<<std::endl;
                std::vector<uint8_t> temp_meta(meta_size);
                meta.resize(meta_size);
                read(meta.data(),meta_size,c);
               // std::cout<<"dwad2"<<std::endl;
                //std::copy(temp_meta.begin(),temp_meta.end(),meta.begin());
                //std::cout<<"dwad3"<<std::endl;
            }
        }

        void print() {
            printf("CmprAlgo = %s\n", enum2Str((ALGO) cmprAlgo));
        }

        char N;
        std::vector<size_t> dims;
        size_t num;
        uint8_t cmprAlgo = ALGO_INTERP_LORENZO;
        uint8_t errorBoundMode = EB_ABS;
        uint8_t tuningTarget=TUNING_TARGET_RD;
        double absErrorBound;
        double relErrorBound=-1.0;
        double psnrErrorBound;
        double l2normErrorBound;
        //double prewave_absErrorBound;
        double rng=-1;
        double alpha=-1;
        double beta=-1;
        double autoTuningRate=0.0;
        double predictorTuningRate=0.0;
        bool lorenzo = true;
        bool lorenzo2 = false;
        bool regression = true;
        bool regression2 = false;
        bool openmp = false;
        uint8_t lossless = 1; // 0-> skip lossless(use lossless_bypass); 1-> zstd
        uint8_t encoder = 1;// 0-> skip encoder; 1->HuffmanEncoder; 2->ArithmeticEncoder
        uint8_t interpAlgo = INTERP_ALGO_CUBIC;
    
        uint8_t interpDirection = 0;
        int levelwisePredictionSelection=0;
        std::vector <uint8_t> interpAlgo_list;
        std::vector <uint8_t> interpDirection_list;
        
        size_t maxStep=0;
        int interpBlockSize = 32;
        int quantbinCnt = 65536;
        int blockSize;
        //int exhaustiveTuning=0;
        int testLorenzo=0;
        std::vector<int> quant_bins;
        //double pred_square_error;
        double decomp_square_error;
        std::vector<size_t> quant_bin_counts;
        int sampleBlockSize=0;
        bool blockwiseTuning=0;
        int stride; //not used now
        int pred_dim; // not used now
        //bool linearReduce=0;
        //bool train=0;
        //bool useCoeff=0;
        //int regSampleStep=6;
        //int multiDimInterp=0;
        int profiling=0;//since there may be multiple ways of profiling set it to int
        int SSIMBlockSize=8;
        int fixBlockSize=0;
        int blockwiseSampleBlockSize=0;
        //std::vector<double> lorenzo1_coeffs;
        //std::vector<double> lorenzo2_coeffs;
        bool verbose=1;
        bool QoZ=0;
        //bool pdTuningRealComp=0;
        //int pdTuningAbConf=0;
        //double pdAlpha=-1;
        //double pdBeta=-1;
        //bool lastPdTuning=0;
        int abList=0;
        int crossBlock=0;
        int sampleBlockSampleBlockSize=0;
        bool peTracking=0;
        int wavelet=0; //may have different wavelets
       //bool external_wave=0;
        std::vector<size_t> coeffs_dims;
        size_t coeffs_num=0;
        double wavelet_rel_coeff = 0.75;
        size_t firstSize;
        int coeffTracking=0;//0 no. 1: output coeff. 2: print stats of coeff 3: both
        int pid=0;

        int offsetPredictor=0;//0:zeropredictor 1: 1D lorenzo 2: MD lorenzo 3:1D interp 4: MD interp
        //int transformation = 0; //0: no trans; 1: sigmoid 2: tanh
        std::vector<float> predictionErrors;//for debug, to delete in final version.
        std::vector<uint8_t> interp_ops;//for debug, to delete in final version.
        std::vector<uint8_t> interp_dirs;//for debug, to delete in final version.
        int trimToZero=0;//1: trim only when quantizing;2: also trim before compression.
        double preTrim=0.0;//trim small numbers to zero before compression.
        int blockOrder = 0;//order of blocks.
        int waveletTest = 1;
        double waveletTuningRate = 0.0;
        int waveletAutoTuning = 0;
        double waveletBrFix = 1.0;
        double waveletMseFix = 1.0;
        double waveletBrFix2 = 1.0;
        double waveletMseFix2 = 1.0;
        double lorenzoBrFix = 1.0;
        int var_first=0;
        size_t profStride=0;
        int sperr=-1;
        //double sperr_eb_coeff = 1.5;
        int waveAutoFix=1;
        int conditioning=0;
        size_t meta_size=0;
        sperr::vec8_type meta;
        std::vector<sperr::vec8_type> block_metas;
        int fixWave=-1;
        bool sperrWithoutWave=false;
        bool pyBind=true;
        std::string metadata;
        bool pybind_activated=false;
        bool FZ=false;
        //bool profilingFix=true;//only for test

       // double anchorThreshold=0.0;
       // size_t minAnchorLevel=3;


        

    };


}

#endif //SZ_CONFIG_HPP

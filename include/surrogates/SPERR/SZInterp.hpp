#ifndef SPERR_SZINTERP_HPP
#define SPERR_SZINTERP_HPP

#include <SPERR/utils/Config.hpp>
#include <SPERR/sperr/SZInterpolationCompressor.hpp>

#include <SPERR/sperr/SZBlockInterpolationCompressor.hpp>

#include <QoZ/preprocessor/Wavelet.hpp>

#include <QoZ/quantizer/IntegerQuantizer.hpp>
#include <QoZ/lossless/Lossless_zstd.hpp>
#include <SPERR/utils/Iterator.hpp>
#include <SPERR/Sample.hpp>
#include <SPERR/utils/Transform.hpp>
#include <SPERR/utils/Statistic.hpp>
#include <SPERR/utils/Extraction.hpp>
#include <SPERR/utils/QuantOptimization.hpp>
#include <SPERR/utils/Metrics.hpp>
#include <SPERR/sperr/ZeroPredictor.hpp>
#include <SPERR/utils/CoeffRegression.hpp>
#include <SPERR/utils/ExtractRegData.hpp>
#include <QoZ/api/impl/SZLorenzoReg.hpp>

#include <SPERR/sperr/SPERR3D_OMP_C.h>
#include <SPERR/sperr/SPERR3D_OMP_D.h>

#include <cmath>
#include <memory>
#include <limits>
#include <cstring>
#include <cstdlib>

#include "qcat_ssim.h"

template<class T, QoZ::uint N>
bool use_sperr(const QoZ::Config & conf){
    return ( (conf.wavelet>0 or conf.sperrWithoutWave) and conf.sperr>=conf.wavelet and N==3);
}

template<class T, QoZ::uint N>
auto pre_Condition(const QoZ::Config &conf,T * data){//conditioner not updated to the newest version of SPERR.
    std::vector<double> buf(data,data+conf.num);//maybe not efficient
    sperr::Conditioner conditioner;
    sperr::dims_type temp_dims={0,0,0};//temp. Fix for later introducing custom filter.
    /*
    if(conf.conditioning==2){
        std::array<bool, 4> b4{true,true,false,false};
        conditioner.toggle_all_settings(b4);
    }
    */

    auto condi_meta = conditioner.condition(buf,temp_dims);
    //if(rtn!=sperr::RTNType::Good)
        //std::cout<<"bad cond"<<std::endl;
    for(size_t i=0;i<conf.num;i++)
        data[i]=buf[i];
    //memcpy(data,buf.data(),conf.num*sizeof(T));//maybe not efficient
    return condi_meta;
}

template<class T, QoZ::uint N>
auto post_Condition(T * data,const size_t &num,const sperr::vec8_type& meta){
    std::vector<double> buf(data,data+num);
   
    sperr::dims_type temp_dims={0,0,0};//temp. Fix for later introducing custom filter.
    sperr::Conditioner conditioner;
    auto rtn = conditioner.inverse_condition(buf,temp_dims,meta);
    for(size_t i=0;i<num;i++)
        data[i]=buf[i];
    //memcpy(data,buf.data(),num*sizeof(T));//maybe not efficient
    return rtn;
}

template<class T, QoZ::uint N> 
char *SPERR_Compress(QoZ::Config &conf, T *data, size_t &outSize){
    assert(N == 1 or N==2 or N==3);//need to complete 2D support later.
        
    SPERR3D_OMP_C compressor;
    compressor.set_num_threads(1);
    compressor.set_eb_coeff(conf.wavelet_rel_coeff);
    if(conf.wavelet!=1)
        compressor.set_skip_wave(true);
    auto rtn = sperr::RTNType::Good;
      
    auto chunks = std::vector<size_t>{1024,1024,1024};//ori 256^3, to tell the truth this is not large enough for scale but I just keep it, maybe set it large later.
    if(N==3)
        rtn = compressor.copy_data(reinterpret_cast<const float*>(data), conf.num,
                                {conf.dims[2], conf.dims[1], conf.dims[0]}, {chunks[0], chunks[1], chunks[2]});
    else if(N==2)
        rtn = compressor.copy_data(reinterpret_cast<const float*>(data), conf.num,
                                {conf.dims[1], conf.dims[0], 1}, {chunks[0], chunks[1], chunks[2]});//temp 2D support. not sure if works well.
    else
        rtn = compressor.copy_data(reinterpret_cast<const float*>(data), conf.num,
                                {conf.dims[0], 1, 1}, {chunks[0], chunks[1], chunks[2]});

    compressor.set_target_pwe(conf.absErrorBound);
    rtn = compressor.compress();
    auto stream = compressor.get_encoded_bitstream();
        
    char * outData=new char[stream.size()+conf.size_est()];
    outSize=stream.size();
    memcpy(outData,stream.data(),stream.size());//maybe not efficient
    stream.clear();
    stream.shrink_to_fit();
    return outData;

}
template<class T, QoZ::uint N> 
void SPERR_Decompress(char *cmpData, size_t cmpSize, T *decData){
    
    std::vector<uint8_t> in_stream(cmpData,cmpData+cmpSize);
    SPERR3D_OMP_D decompressor;
  
    decompressor.set_num_threads(1);
    if (decompressor.use_bitstream(in_stream.data(), in_stream.size()) != sperr::RTNType::Good) {
        std::cerr << "Read compressed file error: "<< std::endl;
        return;
    }

    if (decompressor.decompress(in_stream.data()) != sperr::RTNType::Good) {
        std::cerr << "Decompression failed!" << std::endl;
        return ;
    }
   
    in_stream.clear();
    in_stream.shrink_to_fit();
    const auto vol = decompressor.get_data<float>();
    memcpy(decData,vol.data(),sizeof(T)*vol.size());//maybe not efficient
    return;
}


template<class T, QoZ::uint N>
char * outlier_compress(QoZ::Config &conf,T *data,size_t &outSize){

    char * outlier_compress_output;
    if (conf.offsetPredictor ==0){
        auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
        auto sz = QoZ::make_sz_general_compressor<T, 1>(QoZ::make_sz_general_frontend<T, 1>(conf, QoZ::ZeroPredictor<T, 1>(), quantizer), QoZ::HuffmanEncoder<int>(),
                                                                       QoZ::Lossless_zstd());  
        outlier_compress_output =  (char *)sz->compress(conf,data,outSize);
        delete sz;
    }
    else if (conf.offsetPredictor ==1){
        conf.lorenzo = true;
        conf.lorenzo2 = true;
        conf.regression = false;
        conf.regression2 = false;
        conf.openmp = false;
        conf.blockSize = 16;//original 5
        conf.quantbinCnt = 65536 * 2;

        auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
        auto sz = make_lorenzo_regression_compressor<T, 1>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        outlier_compress_output =  (char *)sz->compress(conf,data,outSize);
        delete sz;
    }
    else if (conf.offsetPredictor == 2){
        conf.setDims(conf.dims.begin(),conf.dims.end());
        conf.lorenzo = true;
        conf.lorenzo2 = true;
        conf.regression = false;
        conf.regression2 = false;
        conf.openmp = false;
        conf.blockSize = 5;
        conf.quantbinCnt = 65536 * 2;
        auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        outlier_compress_output =  (char *)sz->compress(conf,data,outSize);
        delete sz;
    }

    else if (conf.offsetPredictor == 3){
        conf.interpAlgo=QoZ::INTERP_ALGO_CUBIC;
        conf.interpDirection=0;
        auto sz = QoZ::SZInterpolationCompressor<T, 1, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
        QoZ::LinearQuantizer<T>(conf.absErrorBound),
        QoZ::HuffmanEncoder<int>(),
        QoZ::Lossless_zstd());
        outlier_compress_output =  (char *)sz.compress(conf,data,outSize);
        
    }

    else if (conf.offsetPredictor == 4){
            
        conf.setDims(conf.dims.begin(),conf.dims.end());
        conf.interpAlgo=QoZ::INTERP_ALGO_CUBIC;
        conf.interpDirection=0;
        auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
        QoZ::LinearQuantizer<T>(conf.absErrorBound),
        QoZ::HuffmanEncoder<int>(),
        QoZ::Lossless_zstd());
            
       
        outlier_compress_output =  (char *)sz.compress(conf,data,outSize);
        
    }
    return outlier_compress_output;

}

template<class T, QoZ::uint N>
void outlier_decompress(QoZ::Config &conf,char *cmprData,size_t outSize,T*decData){
    if (conf.offsetPredictor ==0){
        auto sz = QoZ::make_sz_general_compressor<T, 1>(QoZ::make_sz_general_frontend<T, 1>(conf, QoZ::ZeroPredictor<T, 1>(), QoZ::LinearQuantizer<T>()), QoZ::HuffmanEncoder<int>(),
                                                                       QoZ::Lossless_zstd());

        sz->decompress((QoZ::uchar *)cmprData,outSize,decData);
       
        delete sz;
    }

    else if (conf.offsetPredictor ==1){
        conf.lorenzo = true;
        conf.lorenzo2 = true;
        conf.regression = false;
        conf.regression2 = false;
        conf.openmp = false;
        conf.blockSize = 16;//original 5
        conf.quantbinCnt = 65536 * 2;

        auto sz = make_lorenzo_regression_compressor<T, 1>(conf, QoZ::LinearQuantizer<T>(), QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
                  
        sz->decompress((QoZ::uchar *)cmprData,outSize,decData);
        delete sz;
    }
    else if (conf.offsetPredictor == 2){
        conf.setDims(conf.dims.begin(),conf.dims.end());
        conf.lorenzo = true;
        conf.lorenzo2 = true;
        conf.regression = false;
        conf.regression2 = false;
        conf.openmp = false;
        conf.blockSize = 5;
        conf.quantbinCnt = 65536 * 2;

        auto sz = make_lorenzo_regression_compressor<T, N>(conf, QoZ::LinearQuantizer<T>(), QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());       
        sz->decompress((QoZ::uchar *)cmprData,outSize,decData);
        delete sz;
    }

    else if (conf.offsetPredictor == 3){
        conf.interpAlgo=QoZ::INTERP_ALGO_CUBIC;
        conf.interpDirection=0;

           
        auto sz = QoZ::SZInterpolationCompressor<T, 1, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
        QoZ::LinearQuantizer<T>(),
        QoZ::HuffmanEncoder<int>(),
        QoZ::Lossless_zstd());      
        sz.decompress((QoZ::uchar *)cmprData,outSize,decData);
        
    }

    else if (conf.offsetPredictor == 4){
            
        conf.setDims(conf.dims.begin(),conf.dims.end());
        conf.interpAlgo=QoZ::INTERP_ALGO_CUBIC;
        conf.interpDirection=0;          
        auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
        QoZ::LinearQuantizer<T>(),
        QoZ::HuffmanEncoder<int>(),
        QoZ::Lossless_zstd());      
        sz.decompress((QoZ::uchar *)cmprData,outSize,decData);
        
    }
    
}

template<class T, QoZ::uint N>
char *SZ_compress_Interp(QoZ::Config &conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interpAlgo << std::endl
//              << "Direction          = " << direction << std::endl
//              << "SZ block size      = " << blockSize << std::endl
//              << "Interp block size  = " << interpBlockSize << std::endl;

    assert(N == conf.N);
    assert(conf.cmprAlgo == QoZ::ALGO_INTERP);
    QoZ::calAbsErrorBound(conf, data);

    //conf.print();
    
    auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
            QoZ::LinearQuantizer<T>(conf.absErrorBound),
            QoZ::HuffmanEncoder<int>(),
            QoZ::Lossless_zstd());

   
    //QoZ::Timer timer;

    //timer.start();
    char *cmpData = (char *) sz.compress(conf, data, outSize);
     //double incall_time = timer.stop();
    //std::cout << "incall time = " << incall_time << "s" << std::endl;
    return cmpData;
}

template<class T, QoZ::uint N>
void SZ_decompress_Interp(QoZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    /*

    assert(conf.cmprAlgo == QoZ::ALGO_INTERP);
    QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;
    if (conf.wavelet==0 and !use_sperr<T,N>(conf)){
        
        auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                QoZ::LinearQuantizer<T>(),
                QoZ::HuffmanEncoder<int>(),
                QoZ::Lossless_zstd());
        if (!conf.blockwiseTuning)
            sz.decompress(cmpDataPos, cmpSize, decData);
        else{
            sz.decompress_block(cmpDataPos, cmpSize, decData);
        }
    }
    
    else{


        if(use_sperr<T,N>(conf) and conf.wavelet<=1){
            std::vector<uint8_t> in_stream(cmpData,cmpData+cmpSize);
            SPERR3D_OMP_D decompressor;
            decompressor.set_num_threads(1);
            if (decompressor.use_bitstream(in_stream.data(), in_stream.size()) != sperr::RTNType::Good) {
                std::cerr << "Read compressed file error: "<< std::endl;
                return;
            }

            if (decompressor.decompress(in_stream.data()) != sperr::RTNType::Good) {
                std::cerr << "Decompression failed!" << std::endl;
                return ;
            }
            in_stream.clear();
            in_stream.shrink_to_fit();
            const auto vol = decompressor.get_data<float>();
            memcpy(decData,vol.data(),sizeof(T)*conf.num);//maybe not efficient
            
            return;




        }
      
        size_t first =conf.firstSize;
        size_t second=cmpSize-conf.firstSize;    

        if(use_sperr<T,N>(conf))
            SPERR_Decompress<T,N>((char*)cmpDataPos, first,decData);
        else{
            auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                    QoZ::LinearQuantizer<T>(),
                    QoZ::HuffmanEncoder<int>(),
                    QoZ::Lossless_zstd());

        
            if (!conf.blockwiseTuning)
                sz.decompress(cmpDataPos, first, decData);
            else{
               
                sz.decompress_block(cmpDataPos, first, decData);
            }
        }
      
        //QoZ::writefile<T>("waved.qoz.dec.sigmo", decData, conf.num);
     

         //QoZ::writefile<T>("waved.qoz.dec.logit", decData, conf.num);

        if(conf.wavelet>1){
            T* newDecData;
            if(conf.pyBind){
              
                
                std::vector<size_t> ori_dims=conf.dims;
                size_t ori_num=conf.num;
                conf.dims=conf.coeffs_dims;
                conf.num=conf.coeffs_num; 
                newDecData= QoZ::pybind_wavelet_postprocessing<T,N>(conf,decData,conf.metadata,conf.wavelet, false,ori_dims);
                conf.dims=ori_dims;
                conf.num=ori_num;
               
                
                
            }
            else
                newDecData= QoZ::external_wavelet_postprocessing<T,N>(decData, conf.coeffs_dims, conf.coeffs_num, conf.wavelet, conf.pid, false,conf.dims);

            
          
            delete []decData;
            decData = new T [conf.num];
            memcpy(decData,newDecData,sizeof(T)*conf.num);//maybe not efficient
            delete []newDecData;

        
        }
        
        else{
            QoZ::Wavelet<T,N> wlt;
            wlt.postProcess_cdf97(decData,conf.dims);
        }
       
        if(conf.conditioning and (!use_sperr<T,N>(conf) or conf.wavelet>1)){
            auto rtn=post_Condition<T,N>(decData,conf.num,conf.meta);
                
        }

       
        //QoZ::writefile<T>("waved.qoz.dec.idwt", decData, conf.num);
       
        
        if(second>0){
            T *offsets =new T [conf.num];
            outlier_decompress<T,N>(conf,(char*)(cmpDataPos+first),second,offsets);
        
        
        //QoZ::writefile<T>("waved.qoz.dec.offset", offsets, conf.num); 
            for(size_t i=0;i<conf.num;i++)
                decData[i]+=offsets[i];//maybe not efficient
            delete [] offsets;
        }
        
    }    
    */
}


template<class T, QoZ::uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                  double eb, int interp_op, int direction_op, int block_size) {
    std::vector<T> data1(data, data + num);
    size_t outSize = 0;
    QoZ::Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpAlgo = interp_op;
    conf.interpDirection = direction_op;
    auto sz = QoZ::SZBlockInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
            QoZ::LinearQuantizer<T>(eb),
            QoZ::HuffmanEncoder<int>(),
            QoZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data1.data(), outSize);
    delete[]cmpData;
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}



/*
template<class T, QoZ::uint N>
int compareWavelets(QoZ::Config &conf, std::vector< std::vector<T> > & sampled_blocks){//This is an unfinished API. Not sure whether useful later.
    size_t sampleBlockSize=conf.sampleBlockSize;
    std::vector<size_t> global_dims=conf.dims;
    size_t global_num=conf.num;

    num_sampled_blocks=sampled_blocks.size();
    per_block_ele_num=pow(sampleBlockSize+1,N);
    ele_num=num_sampled_blocks*per_block_ele_num;
    conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
    conf.num=per_block_ele_num;
    std::vector<T> cur_block(per_block_ele_num,0);

    double wave_eb=conf.absErrorBound*conf.wavelet_rel_coeff;

    size_t sig_count=0;

    std::vector<T> gathered_coeffs;
    std::vector<T> gathered_blocks;

    return 0;

}
*/


template<class T, QoZ::uint N>
void sampleBlocks(T *data,std::vector<size_t> &dims, size_t sampleBlockSize,std::vector< std::vector<T> > & sampled_blocks,double sample_rate,int profiling ,std::vector<std::vector<size_t> > &starts,int var_first=0){
    for(int i=0;i<sampled_blocks.size();i++){
                std::vector< T >().swap(sampled_blocks[i]);                
            }
            std::vector< std::vector<T> >().swap(sampled_blocks);
    for(int i=0;i<sampled_blocks.size();i++){
        std::vector< T >().swap(sampled_blocks[i]);                  
    }
    std::vector< std::vector<T> >().swap(sampled_blocks);                               
    size_t totalblock_num=1;
    for(int i=0;i<N;i++){                        
        totalblock_num*=(int)((dims[i]-1)/sampleBlockSize);
    }               
    size_t idx=0,block_idx=0;   
    if(profiling){
        size_t num_filtered_blocks=starts.size();    
        if(var_first==0){  
            size_t sample_stride=(size_t)(num_filtered_blocks/(totalblock_num*sample_rate));
            if(sample_stride<=0)
                sample_stride=1;
            
            for(size_t i=0;i<num_filtered_blocks;i+=sample_stride){
                std::vector<T> s_block;
                QoZ::sample_blocks<T,N>(data, s_block,dims, starts[i],sampleBlockSize+1);
                sampled_blocks.push_back(s_block);
                
            }
            
        }
        else{
            std::vector< std::pair<double,std::vector<size_t> > >block_heap;
            for(size_t i=0;i<num_filtered_blocks;i++){
                double mean,sigma2,range;
                QoZ::blockwise_profiling<T>(data,dims, starts[i],sampleBlockSize+1, mean,sigma2,range);
                block_heap.push_back(std::pair<double,std::vector<size_t> >(sigma2,starts[i]));
                
            }
            std::make_heap(block_heap.begin(),block_heap.end());
          

            size_t sampled_block_num=totalblock_num*sample_rate;
            if(sampled_block_num>num_filtered_blocks)
                sampled_block_num=num_filtered_blocks;
            if(sampled_block_num==0)
                sampled_block_num=1;

            for(size_t i=0;i<sampled_block_num;i++){
                std::vector<T> s_block;
             
                QoZ::sample_blocks<T,N>(data, s_block,dims, block_heap.front().second,sampleBlockSize+1);
              
                sampled_blocks.push_back(s_block);
                std::pop_heap(block_heap.begin(),block_heap.end());
                block_heap.pop_back();
               
            }
        }
    }               
    else{
        if(var_first==0){
            size_t sample_stride=(size_t)(1.0/sample_rate);
            if(sample_stride<=0)
                sample_stride=1;
            if(N==1) {
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){       
                    if (idx%sample_stride==0){
                        std::vector<size_t> starts{x_start};
                        std::vector<T> s_block;
                        QoZ::sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                        sampled_blocks.push_back(s_block);
                    }
                    idx+=1;
                }
            }
            else if (N==2){                        
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                           
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                        if (idx%sample_stride==0){
                            std::vector<size_t> starts{x_start,y_start};
                            std::vector<T> s_block;
                            QoZ::sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                            sampled_blocks.push_back(s_block);
                        }
                        idx+=1;
                    }
                }
            }
            else if (N==3){                  
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                          
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                        for (size_t z_start=0;z_start<dims[2]-sampleBlockSize;z_start+=sampleBlockSize){
                            if (idx%sample_stride==0){
                                std::vector<size_t> starts{x_start,y_start,z_start};
                                std::vector<T> s_block;
                                QoZ::sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                                sampled_blocks.push_back(s_block);
                            }
                            idx+=1;
                        }
                    }
                }
            }
        }
        else{
            std::vector <std::vector<size_t> > blocks_starts;
            if (N==1){  
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){  
                    blocks_starts.push_back(std::vector<size_t>{x_start});
                }

            }
            else if (N==2){  
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                           
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                       
                            blocks_starts.push_back(std::vector<size_t>{x_start,y_start});
                    }
                }

            }
            else if (N==3){           
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                          
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                        for (size_t z_start=0;z_start<dims[2]-sampleBlockSize;z_start+=sampleBlockSize){
                            blocks_starts.push_back(std::vector<size_t>{x_start,y_start,z_start});
                        }
                    }
                }
            

                std::vector< std::pair<double,std::vector<size_t> > >block_heap;
                for(size_t i=0;i<totalblock_num;i++){
                    double mean,sigma2,range;
                    QoZ::blockwise_profiling<T>(data,dims, blocks_starts[i],sampleBlockSize+1, mean,sigma2,range);
                    block_heap.push_back(std::pair<double,std::vector<size_t> >(sigma2,blocks_starts[i]));
                }
                std::make_heap(block_heap.begin(),block_heap.end());
                size_t sampled_block_num=totalblock_num*sample_rate;
                if(sampled_block_num==0)
                    sampled_block_num=1;
                for(size_t i=0;i<sampled_block_num;i++){
                    std::vector<T> s_block;
                    QoZ::sample_blocks<T,N>(data, s_block,dims, block_heap.front().second,sampleBlockSize+1);
                    sampled_blocks.push_back(s_block);
                    std::pop_heap(block_heap.begin(),block_heap.end());
                    block_heap.pop_back();
                }

            }
        }
    }
}


template<class T, QoZ::uint N>
double SSIMTest(const QoZ::Config &conf,const std::vector< std::vector<T> > & sampled_blocks){
    assert(N==1 or N==2 or N==3);
    QoZ::Config testConfig(conf);
    double square_error=0.0;
    T maxValue = sampled_blocks[0][0];
    double bitrate=0.0;
    double metric=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;
    std::vector<T> cur_block(testConfig.num,0);
    size_t idx=0;   
    size_t totalOutSize=0;
    float ssim_sum = 0.0;
    size_t num_ssim_samples = 0;

    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        std::vector<T> cur_block(testConfig.num);
       
        std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());

        char* cmprData = SPERR_Compress<T,N>(testConfig,cur_block.data(),sampleOutSize);
        T* decData = new T[per_block_ele_num];
        SPERR_Decompress<T,N>(cmprData, sampleOutSize, decData);


        double ssim_score;
        // double * min_ssim = new double();
        // double * avg_ssim = new double();
        // double * max_ssim = new double();
        T* blockdata = new T[per_block_ele_num];

        memcpy(blockdata, sampled_blocks[k].data(), sizeof(T)*per_block_ele_num);
        if(N==1){
            ssim_score = SSIM_1d_windowed_float(blockdata, decData, sampleBlockSize, 8, 8);
            if(std::isnan(ssim_score)){
                ssim_score = 0;
            }
            else{
                num_ssim_samples++;
            }
        }
        else if(N == 2){
            // ssim_score = zc_calc_ssim_2d_float(blockdata, decData, sampleBlockSize, sampleBlockSize);
            ssim_score = SSIM_2d_windowed_float(blockdata, decData, sampleBlockSize, sampleBlockSize, 8, 8, 8, 8);
            if(std::isnan(ssim_score)){
                ssim_score = 0;
            }
            else{
                num_ssim_samples++;
            }
        }
        else {
            // zc_calc_ssim_3d_float(blockdata, decData, sampleBlockSize, sampleBlockSize, sampleBlockSize, min_ssim, avg_ssim, max_ssim);
            // ssim_score = (float) *avg_ssim;
            ssim_score = SSIM_3d_windowed_float(blockdata, decData, sampleBlockSize, sampleBlockSize, sampleBlockSize, 8,8,8, 8,8,8);

            // indicates block was all constant, just drop this block from the estimation
            if(std::isnan(ssim_score)){
                ssim_score = 0;
            }
            else{
                num_ssim_samples++;
            }

            // std::cout << "ssim: " << *min_ssim << " " << *avg_ssim << " " << *max_ssim << std::endl;
        }

        ssim_sum += ssim_score;

        // delete min_ssim;
        // delete avg_ssim;
        // delete max_ssim;
        delete[] blockdata;
        
    }
                
    double ssim = (double) ssim_sum / num_ssim_samples;

    size_t num_nan_blocks = num_sampled_blocks - num_ssim_samples;
    std::cout << "NDIM: " << N << std::endl;
    std::cout << "sampleblocksize: " << sampleBlockSize << std::endl;
    std::cout << "Num total blocks: " << num_sampled_blocks << std::endl;
    std::cout << "Num NAN blocks: " << num_nan_blocks << std::endl;

    return ssim;
}

template<class T, QoZ::uint N> 
double estimateSPERRSSIMbasedonErrorBound(double error_bound,T * data, double sample_rate, size_t blocksize,std::vector<size_t> &dims,bool skip_outlier=false){//, int profiling=0,int var_first=0){

    //error_bound: The abs error bound.
    //data: The whole input data array. 
    //sample rate: The rate of sampled data points (e.g. 0.01).
    //blocksize: The sampled block size (e.g. 32).
    //dims: the dimensions of the input data array (fastest last, e.g. for miranda data it is {256,384,384}).
    //skip_outlier: whether skip the outlier correction (dewave and error correction).

    std::vector< std::vector<T> > sampled_blocks;
   
    //size_t num_sampled_blocks;
    //size_t per_block_ele_num;
    //size_t ele_num;

    QoZ::Config conf(1);//maybe a better way exist?
    conf.setDims(dims.begin(),dims.end());
    conf.sperr=1;
    conf.wavelet=1;
    conf.wavelet_rel_coeff=1.5;
    //conf.profiling=profiling;
    //conf.var_first=var_first;
    conf.sampleBlockSize=blocksize;
    conf.cmprAlgo=QoZ::ALGO_INTERP;
    conf.tuningTarget=QoZ::TUNING_TARGET_CR;
    conf.errorBoundMode=QoZ::EB_ABS;
    conf.absErrorBound=error_bound;
    if(skip_outlier)
        conf.wavelet=2;
    //if (conf.rng<0)
     //   conf.rng=QoZ::data_range<T>(data,conf.num);
    size_t totalblock_num=1;  
    for(int i=0;i<N;i++){                      
        totalblock_num*=(size_t)((conf.dims[i]-1)/conf.sampleBlockSize);
    }


    std::vector<std::vector<size_t> >starts;
    if(conf.profiling){      
        conf.profStride=conf.sampleBlockSize/4;
        if(N==1){
            QoZ::profiling_block_1d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
        else if(N==2){
            QoZ::profiling_block_2d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
        else if (N==3){
            QoZ::profiling_block_3d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
       
    }
  

    size_t num_filtered_blocks=starts.size();

    sampleBlocks<T,N>(data,conf.dims,conf.sampleBlockSize,sampled_blocks,sample_rate,conf.profiling,starts,conf.var_first);
   
           
    //num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(blocksize+1,N);
   // ele_num=num_sampled_blocks*per_block_ele_num;

    std::vector<size_t> blockdims(N,blocksize+1);
    conf.setDims(blockdims.begin(),blockdims.end());
   


    double est_ssim = SSIMTest<T,N>(conf, sampled_blocks);

    return est_ssim;


    /*
    conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
    conf.num=per_block_ele_num;
    std::vector<T> cur_block(per_block_ele_num,0);
    */

}

template<class T, QoZ::uint N>
double PSNRTest(const QoZ::Config &conf,const std::vector< std::vector<T> > & sampled_blocks, T value_range){
    QoZ::Config testConfig(conf);
    double square_error=0.0;
    T maxValue = sampled_blocks[0][0];
    T minValue = sampled_blocks[0][0];
    double bitrate=0.0;
    double metric=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;
    std::vector<T> cur_block(testConfig.num,0);
    size_t idx=0;   
    size_t totalOutSize=0;

    double estPSNR = 0.0;
    size_t num_psnr_samples = 0;
    
    if(num_sampled_blocks == 0){
        std::cerr << "Warning: Num sampled blocks is zero! Consider using a smaller blocksize."<< std::endl;
    }
                           
    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        std::vector<T> cur_block(testConfig.num);
       
        std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());

        char* cmprData = SPERR_Compress<T,N>(testConfig,cur_block.data(),sampleOutSize);
        T* decData = new T[per_block_ele_num];
        SPERR_Decompress<T,N>(cmprData, sampleOutSize, decData);

        // square_error = 0.0;
        maxValue = sampled_blocks[k][0];
        for(int j = 0; j < per_block_ele_num; j++){
            square_error += pow(decData[j] - sampled_blocks[k][j], 2);
            if(sampled_blocks[k][j] > maxValue){
                maxValue = sampled_blocks[k][j];
            }
            if(sampled_blocks[k][j] < minValue){
                minValue = sampled_blocks[k][j];
            }
        }

        // double eps = 1e-16;
        // double mse = square_error / per_block_ele_num;
        // // T value_range = maxValue - minValue;
        // double psnr = -20.0*log10((sqrt(mse) / value_range) + eps);

        

        // std::cout << square_error << " " << mse << " " << psnr << std::endl;

        // if(std::isnan(psnr)){
        //     // no block error, simulate high psnr
        //     // std::cout << "max " << maxValue << std::endl;
        //     psnr = 320.0;
        // }
        
        // estPSNR += psnr;
        // num_psnr_samples++;
        

        delete[] decData;
        delete[] cmprData;
        
                
    }
                
    double eps = 1e-16;
    double mse = square_error / ele_num;
	estPSNR = -20.0*log10((sqrt(mse) / value_range) + eps);   
    // std::cout << estPSNR << "\n";  
    // estPSNR /= num_psnr_samples;
    // std::cout << estPSNR << std::endl;  

    // std::cout << "Num blocks dropped / total: " << (num_sampled_blocks - num_psnr_samples) << " / " << num_sampled_blocks << std::endl;

    return estPSNR;
}

template<class T, QoZ::uint N> 
double estimateSPERRPSNRbasedonErrorBound(double error_bound,T * data, double sample_rate, size_t blocksize,std::vector<size_t> &dims,bool skip_outlier=false){//, int profiling=0,int var_first=0){

    //error_bound: The abs error bound.
    //data: The whole input data array. 
    //sample rate: The rate of sampled data points (e.g. 0.01).
    //blocksize: The sampled block size (e.g. 32).
    //dims: the dimensions of the input data array (fastest last, e.g. for miranda data it is {256,384,384}).
    //skip_outlier: whether skip the outlier correction (dewave and error correction).

    std::vector< std::vector<T> > sampled_blocks;
   
    //size_t num_sampled_blocks;
    //size_t per_block_ele_num;
    //size_t ele_num;

    QoZ::Config conf(1);//maybe a better way exist?
    conf.setDims(dims.begin(),dims.end());
    conf.sperr=1;
    conf.wavelet=1;
    conf.wavelet_rel_coeff=1.5;
    //conf.profiling=profiling;
    //conf.var_first=var_first;
    conf.sampleBlockSize=blocksize;
    conf.cmprAlgo=QoZ::ALGO_INTERP;
    conf.tuningTarget=QoZ::TUNING_TARGET_CR;
    conf.errorBoundMode=QoZ::EB_ABS;
    conf.absErrorBound=error_bound;
    if(skip_outlier)
        conf.wavelet=2;
    //if (conf.rng<0)
     //   conf.rng=QoZ::data_range<T>(data,conf.num);
    size_t totalblock_num=1;  
    for(int i=0;i<N;i++){                      
        totalblock_num*=(size_t)((conf.dims[i]-1)/conf.sampleBlockSize);
    }

    T max = data[0];
    T min = data[0];
    for(int i = 0; i < conf.num; i++){
        if(data[i]>max) max = data[i];
        if(data[i]<min) min=data[i];
    }


    std::vector<std::vector<size_t> >starts;
    if(conf.profiling){      
        conf.profStride=conf.sampleBlockSize/4;
        if(N==1){
            QoZ::profiling_block_1d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
        else if(N==2){
            QoZ::profiling_block_2d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
        else if (N==3){
            QoZ::profiling_block_3d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
       
    }
  

    size_t num_filtered_blocks=starts.size();

    sampleBlocks<T,N>(data,conf.dims,conf.sampleBlockSize,sampled_blocks,sample_rate,conf.profiling,starts,conf.var_first);
   
           
    //num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(blocksize+1,N);
   // ele_num=num_sampled_blocks*per_block_ele_num;

    std::vector<size_t> blockdims(N,blocksize+1);
    conf.setDims(blockdims.begin(),blockdims.end());
   

    T value_range = max - min;
    double est_psnr = PSNRTest<T,N>(conf, sampled_blocks, value_range);

    return est_psnr;


    /*
    conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
    conf.num=per_block_ele_num;
    std::vector<T> cur_block(per_block_ele_num,0);
    */

}

template<class T, QoZ::uint N>
std::pair<double,double> CompressTest(const QoZ::Config &conf,const std::vector< std::vector<T> > & sampled_blocks){
    QoZ::Config testConfig(conf);
    double square_error=0.0;
    double bitrate=0.0;
    double metric=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;
    std::vector<T> cur_block(testConfig.num,0);
    size_t idx=0;   
    size_t totalOutSize=0;
                           
    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        std::vector<T> cur_block(testConfig.num);
       
        std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());

        auto cmprData=SPERR_Compress<T,N>(testConfig,cur_block.data(),sampleOutSize);
        
        totalOutSize+=sampleOutSize;

        double curr_bitrate=8*double(sampleOutSize)/per_block_ele_num;
        if(testConfig.wavelet==1){
            curr_bitrate*=testConfig.waveletBrFix;
        } 
        else if(testConfig.wavelet>1){
            curr_bitrate*=testConfig.waveletBrFix2;
        }  

        bitrate += curr_bitrate;
                
    }

    // take average bitrate over blocks
    bitrate = bitrate / num_sampled_blocks;
                
    // bitrate=8*double(totalOutSize)/ele_num;
    // if(testConfig.wavelet==1){
    //     bitrate*=testConfig.waveletBrFix;
    // } 
    // else if(testConfig.wavelet>1){
    //     bitrate*=testConfig.waveletBrFix2;
    // }       

    //std::cout<<bitrate<<" "<<metric<<std::endl;
    return std::pair(bitrate,metric);
}



template<class T, QoZ::uint N> 
double estimateSPERRCRbasedonErrorBound(double error_bound,T * data, double sample_rate, size_t blocksize,std::vector<size_t> &dims,bool skip_outlier=false){//, int profiling=0,int var_first=0){

    //error_bound: The abs error bound.
    //data: The whole input data array. 
    //sample rate: The rate of sampled data points (e.g. 0.01).
    //blocksize: The sampled block size (e.g. 32).
    //dims: the dimensions of the input data array (fastest last, e.g. for miranda data it is {256,384,384}).
    //skip_outlier: whether skip the outlier correction (dewave and error correction).

    std::vector< std::vector<T> > sampled_blocks;
   
    //size_t num_sampled_blocks;
    //size_t per_block_ele_num;
    //size_t ele_num;

    QoZ::Config conf(1);//maybe a better way exist?
    conf.setDims(dims.begin(),dims.end());
    conf.sperr=1;
    conf.wavelet=1;
    conf.wavelet_rel_coeff=1.5;
    //conf.profiling=profiling;
    //conf.var_first=var_first;
    conf.sampleBlockSize=blocksize;
    conf.cmprAlgo=QoZ::ALGO_INTERP;
    conf.tuningTarget=QoZ::TUNING_TARGET_CR;
    conf.errorBoundMode=QoZ::EB_ABS;
    conf.absErrorBound=error_bound;
    if(skip_outlier)
        conf.wavelet=2;
    //if (conf.rng<0)
     //   conf.rng=QoZ::data_range<T>(data,conf.num);
    size_t totalblock_num=1;  
    for(int i=0;i<N;i++){                      
        totalblock_num*=(size_t)((conf.dims[i]-1)/conf.sampleBlockSize);
    }


    std::vector<std::vector<size_t> >starts;
    if(conf.profiling){      
        conf.profStride=conf.sampleBlockSize/4;
        if(N==1){
            QoZ::profiling_block_1d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
        else if(N==2){
            QoZ::profiling_block_2d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
        else if (N==3){
            QoZ::profiling_block_3d<T,N>(data,conf.dims,starts,blocksize,conf.absErrorBound,conf.profStride);
        }
       
    }
  

    size_t num_filtered_blocks=starts.size();

    sampleBlocks<T,N>(data,conf.dims,conf.sampleBlockSize,sampled_blocks,sample_rate,conf.profiling,starts,conf.var_first);
   
           
    //num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(blocksize+1,N);
   // ele_num=num_sampled_blocks*per_block_ele_num;

    std::vector<size_t> blockdims(N,blocksize+1);
    conf.setDims(blockdims.begin(),blockdims.end());
   


    std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks);

    double cur_ratio=sizeof(T)*8.0/results.first;
  
    

    return cur_ratio;


    /*
    conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
    conf.num=per_block_ele_num;
    std::vector<T> cur_block(per_block_ele_num,0);
    */

}





#endif
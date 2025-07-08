#ifndef SZ3_WAVELET_HPP
#define SZ3_WAVELET_HPP


#include "QoZ/preprocessor/PreProcessor.hpp"
#include "QoZ/sperr/CDF97.h"
#ifdef ENABLE_GSL 
#include <gsl/gsl_wavelet.h>
#endif
#include "QoZ/utils/FileUtil.hpp"
//#include <pybind11/embed.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
//namespace py = pybind11;
namespace QoZ {


    template<class T, QoZ::uint N>
    T * external_wavelet_preprocessing(T *data, const std::vector<size_t> &dims, size_t num, int wave_type=2, size_t pid=0, bool inplace=true, std::vector<size_t> &coeffs_size=std::vector<size_t>())
    {
        std::string input_filename = std::to_string(pid) + "_external_wave_temp_input.tmp";
        //std::cout<<num<<std::endl;
        std::string del_command="rm -f "+ input_filename;
        system(del_command.c_str());
        QoZ::writefile<T>(input_filename.c_str(), data, num);


        std::string wavetype;
        /*Interp
        if (wave_type==2)
            wavetype="sym16";//rtms
        else if(wave_type==3)
            wavetype="bior3.1";//miranda scale
        else if(wave_type==4)
            wavetype="bior4.4";//nyx (hurricane)
        else if(wave_type==5)
            wavetype="coif6";//qmcpack
        else
            wavetype="bior6.8";//hurricane
        */

        
        if (wave_type==2)
            wavetype="sym13";//rtms,hurricane,nyx,miranda
        else if(wave_type==3)
            wavetype="bior3.3";//qmcpack
        else if(wave_type==4)
            wavetype="bior4.4";//scale
        

     

        //std::cout<<wavetype<<std::endl;

       
        std::string command = "python coeff_dwt.py " + input_filename + " " + wavetype + " " + std::to_string(pid);
        for (int i = N - 1; i >= 0; i--)
        {
            command += " " + std::to_string(dims[i]);
        }

        system(command.c_str());
        system(del_command.c_str());

        

        std::string coeffs_filename = std::to_string(pid) + "_external_wave_coeffs.tmp";
        del_command="rm -f "+ coeffs_filename;

        if (inplace)
        {
            QoZ::readfile<T>(coeffs_filename.c_str(), num, data);
            system(del_command.c_str());
            return data;
        }
        else
        {
            coeffs_size.resize(N);
            std::string size_filename = std::to_string(pid) + "_external_coeffs_size.tmp";
            QoZ::readfile<size_t>(size_filename.c_str(), N, coeffs_size.data());


            //for (int i = 0; i <N; i++)
            //{
                //std::cout<<coeffs_size[i]<<std::endl;
           // }

            size_t coeffs_num = 1;
            for (size_t i = 0; i < N; i++)
                coeffs_num *= coeffs_size[i];
            //std::cout<<coeffs_num<<std::endl;

            T *coeffData = new T[coeffs_num];
            QoZ::readfile<T>(coeffs_filename.c_str(), coeffs_num, coeffData);
            system(del_command.c_str());
            del_command="rm -f "+ size_filename;
            system(del_command.c_str());
            return coeffData;
        }
    }

    template<class T, QoZ::uint N>
    T * external_wavelet_postprocessing(T *data, const std::vector<size_t> &dims, size_t num, int wave_type=2, size_t pid=0, bool inplace=true,const std::vector<size_t> &output_dims=std::vector<size_t>())
    {
        
            
        std::string input_filename = std::to_string(pid) + "_external_wave_coeff_input.tmp";
        std::string del_command="rm -f "+ input_filename;
        system(del_command.c_str());
        
        QoZ::writefile<T>(input_filename.c_str(), data, num);
        std::string command = "python coeff_idwt.py " + input_filename;
                
        system(command.c_str());
       
        system(del_command.c_str());
        std::string output_filename = std::to_string(pid) + "_external_deccoeff_idwt.tmp";
        del_command="rm -f "+ output_filename;
        
      
        if (inplace)
        {
            QoZ::readfile<T>(output_filename.c_str(), num, data);
            
            system(del_command.c_str());
            return data;
        }
        else
        {
            size_t outnum=1;
            for (size_t i = 0; i < N; i++)
                outnum *= output_dims[i];

            T *outData = new T[outnum];
            QoZ::readfile<T>(output_filename.c_str(), outnum, outData);
            system(del_command.c_str());
            
            return outData;
        }
    }


    
   
    template<class T, uint N>

    class Wavelet : public concepts::PreprocessorInterface<T, N> {
    public:
        

        #ifdef ENABLE_GSL
        void preProcess(T *data, size_t n) {




            
            size_t m = n - 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m++;

            std::vector<double> dwtdata(m, 0);
            gsl_wavelet *w;
            gsl_wavelet_workspace *work;

            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
            work = gsl_wavelet_workspace_alloc(m);

            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            int status = gsl_wavelet_transform_forward(w, dwtdata.data(), 1, m, work);

            if (status != GSL_SUCCESS) {
                printf("Error: wavelets transform failed.\n");
                exit(0);
            }

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }

            gsl_wavelet_free(w);
            gsl_wavelet_workspace_free(work);
            

        }
        
        void postProcess(T *data, size_t n) {
            size_t m = n - 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m++;

            std::vector<double> dwtdata(m, 0);
            gsl_wavelet *w;
            gsl_wavelet_workspace *work;

            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
            work = gsl_wavelet_workspace_alloc(m);

            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            int status = gsl_wavelet_transform_inverse(w, dwtdata.data(), 1, m, work);

            if (status != GSL_SUCCESS) {
                printf("Error: wavelets transform failed.\n");
                exit(0);
            }

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }

            gsl_wavelet_free(w);
            gsl_wavelet_workspace_free(work);

        }
        #endif
        void preProcess_cdf97(T *data, std::vector<size_t> dims) {
            size_t n=1;
            std::array<size_t,3> m_dims=std::array<size_t,3>{1,1,1};
            for (size_t i=0;i<N;i++){
                n*=dims[i];
                m_dims[N-1-i]=dims[i];
            }

            std::vector<double> dwtdata(n, 0);
            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            sperr::CDF97 m_cdf;

            m_cdf.take_data(std::move(dwtdata), m_dims);
            /*
            auto xforms_xy = num_of_xforms(std::min(m_dims[0], m_dims[1]));
            auto xforms_z = num_of_xforms(m_dims[2]);
            if (xforms_xy == xforms_z)
                m_cdf.dwt3d_dyadic();
            else
                m_cdf.dwt3d_wavelet_packet();
            */
            m_cdf.dwt3d();


            dwtdata=m_cdf.release_data();

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }
            //std::cout<<"pre finished"<<std::endl;
        }


        void postProcess_cdf97(T *data, std::vector<size_t> dims) {
            size_t n=1;
            std::array<size_t,3> m_dims=std::array<size_t,3>{1,1,1};
            for (size_t i=0;i<N;i++){
                n*=dims[i];
                m_dims[N-1-i]=dims[i];
            }
            std::vector<double> dwtdata(n, 0);
            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            sperr::CDF97 m_cdf;

            m_cdf.take_data(std::move(dwtdata), m_dims);
            /*
            auto xforms_xy = num_of_xforms(std::min(m_dims[0], m_dims[1]));
            auto xforms_z = num_of_xforms(m_dims[2]);
            if (xforms_xy == xforms_z)
                m_cdf.idwt3d_dyadic();
            else
                m_cdf.idwt3d_wavelet_packet();
            */
            m_cdf.idwt3d();


            dwtdata=m_cdf.release_data();

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }
            //std::cout<<"post finished"<<std::endl;




        }
        
        
    

       

        
    };
}


#endif //SZ3_WAVELET_HPP

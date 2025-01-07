# Implementing a Compressor Surrogate

ROCCI aims to make implementing a compressor surrogate for autotuning as simple as possible by leveraging `libpressio`, a generic lossy compressor interface. We use `libpressio_opt` to search over configuration settings for each compressor using its corresponding surrogate.

To implement a new compressor surrogate in ROCCI, you only need to do two things: implement your surrogate and register it with ROCCI.

## Implementing a surrogate using `libpressio`

Each compressor surrogate is implemented as a compressor plugin for `libpressio`. This enables us to use `libpressio`'s generic interface to estimate quantities of interest such as compression ratio, psnr, and ssim using eachs surrogate. Please refer to the surrogate implementations in `src/surrogates/` for examples. 

The options and settings your surrogate supports will vary according to the nature of the compressor whose behavior you are estimating and the sophistication of your surrogate model. Your surrogate should perform metric estimation during calls to `compress_impl`. Additionally, your surrogate should report these metrics in `get_metrics_results`. The metrics returned by your compressor will be used by `libpressio_opt` to optimize the parameter configuration of your compressor.

Let's look at the SZx surrogate as an example, in the `compress_impl` function it estimates the appropriate metric based on its configuration.

```c++
int compress_impl(const pressio_data* input,
                struct pressio_data* output) override
{
    double eb = absErrBound;
    if(errBoundMode == REL) {
        eb = relBoundRatio;
        printf("szx_surrogate does not support a relative error bound\n");
        exit(0);
    }
    float* data = (float*) input->data();
    if (metric_setting == CR_METRIC) {
        estCR = szx_estimate_cr_float(data, eb, samplingStride, input->num_elements());
    } else if (metric_setting == PSNR_METRIC) {
        estPSNR = szx_estimate_psnr_float(data, eb, samplingStride, input->num_elements());
    } else if (metric_setting == SSIM_METRIC) {
        estSSIM = szx_estimate_ssim_float(data, eb, samplingStride, blockSize, const_cast<size_t*>(input->dimensions().data()), input->num_elements());
    } else {
        return set_error(1, "invalid metric string, expected one of {cr, psnr, ssim}");
    }

    return 0;

}
```

And then returns these estimated quantities in `get_metrics_results`:

```c++
pressio_options get_metrics_results_impl() const override {
    struct pressio_options options;
    set(options, "szx_surrogate:cr", estCR);
    set(options, "szx_surrogate:psnr", estPSNR);
    set(options, "szx_surrogate:ssim", estSSIM);

    return options;
}
```

These metrics are then used in `roccifastSearchBestSetting` for autotuning.

### Useful resources 
Libpressio is extensively documented here: https://robertu94.github.io/libpressio/

Additional information on implementing a compressor plugin is available here: https://github.com/robertu94/libpressio/blob/master/docs/WritingACompressorPlugin.md

And here: https://github.com/FTHPC/libpressio_tutorial/tree/master/exercises/5_writing_basic_compressors

## Registering your surrogate with ROCCI

Registering a surrogate is simple, just add your surrogate's prefix string to the list in `get_surrogate` in `rocci_utils.cc`.

Make sure to add your surrogate's sources to the root `CMakeLists.txt`.
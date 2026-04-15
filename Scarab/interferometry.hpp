#pragma once
#include "ExternalLibraries/pocketfft_hdronly.h"
#include "batch_acquisition.hpp"
#include "spectrometer.hpp"
#include <chrono>
#include <complex>
#include <format>
#include <fstream>
#include <memory>
#include <span>
#include <vector>

class SpectralInterferometry {
    std::vector<double> wavelengths;
    std::vector<double> frequencies;
    std::vector<double> wavelengths_on_frequency_grid;
    std::vector<double> interference_data;
    std::vector<double> interference_data_interpolated;
    std::vector<double> reference_data_a;
    std::vector<double> reference_data_a_interpolated;
    std::vector<double> reference_data_b;
    std::vector<double> reference_data_b_interpolated;
    std::vector<double> spectral_phase;
    std::vector<double> spectral_phase_mean;
    std::vector<double> spectral_phase_mean_sqr;
    std::vector<double> group_delay;
    size_t phase_count = 0;
    std::vector<std::complex<double>> fft_data;
    std::vector<std::complex<double>> fft_reference_a;
    std::vector<std::complex<double>> fft_reference_b;
    std::vector<double> fft_data_real;
    std::vector<double> fft_reference_a_real;
    std::vector<double> fft_reference_b_real;
    std::vector<double> fft_scale;
    std::vector<double> time_filter;
    std::vector<std::complex<double>> hilbert_time_buffer;
    std::vector<std::complex<double>> hilbert_time_buffer2;
    std::vector<double> hilbert_real_buffer;
    std::vector<double> hilbert_imag_buffer;
    bool is_configured = false;
    bool use_averaging = true;
    double d_f = 2e12;
    double min_f = 100e12;
    double max_f = 2e12 * 512 + 100e12;
    double filter_t0 = 90;
    double filter_width = 50;
    double filter_order = 8;
    size_t num_freqs = 512;
    size_t fft_size = 0;
    void unwrap(std::vector<double> &phase_in_out) {
        for(size_t i = 1; i < num_freqs; i++) {
            double delta = phase_in_out[i] - phase_in_out[i - 1];
            delta -= twoPi<double>() * round(delta / twoPi<double>());
            phase_in_out[i] = phase_in_out[i - 1] + delta;
        }
    }
    void setup_fft() {
        std::unique_ptr<double[]> setup_work_d(new double[num_freqs]);
        std::unique_ptr<std::complex<double>[]> setup_work_c(new std::complex<double>[num_freqs]);
        fft_size = num_freqs;
        hilbert_time_buffer = std::vector<std::complex<double>>(num_freqs / 2 + 1);
        hilbert_time_buffer2 = hilbert_time_buffer;
        hilbert_real_buffer = std::vector<double>(num_freqs);
        hilbert_imag_buffer = hilbert_real_buffer;
        fft_reference_a = hilbert_time_buffer;
        fft_reference_b = hilbert_time_buffer;
        time_filter = std::vector<double>(num_freqs / 2 + 1);
        fft_scale = time_filter;
        fft_reference_a_real = time_filter;
        fft_reference_b_real = time_filter;
        fft_data_real = time_filter;
    }

    void filtered_hilbert(std::vector<double> &in_data,
                          std::vector<double> &out_data_real,
                          std::vector<double> &out_data_imag) {
        if(fft_size != num_freqs) {
            setup_fft();
        }

        pocketfft::r2c(pocketfft::shape_t{in_data.size()},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
                       pocketfft::shape_t{0},
                       pocketfft::FORWARD,
                       in_data.data(),
                       hilbert_time_buffer2.data(),
                       1.0);

        double dt = 1.0e3 / (num_freqs * d_f);
        std::complex<double> ii(0.0, 1.0);
        for(size_t i = 0; i < (num_freqs / 2 + 1); i++) {
            hilbert_time_buffer2[i] *= std::exp(
                -std::pow((static_cast<double>(i) * dt - filter_t0) / filter_width, filter_order) /
                sqrt(2.0));
            hilbert_time_buffer[i] = ii * hilbert_time_buffer2[i];
        }

        pocketfft::c2r(pocketfft::shape_t{out_data_real.size()},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
                       pocketfft::shape_t{0},
                       pocketfft::BACKWARD,
                       hilbert_time_buffer.data(),
                       out_data_real.data(),
                       1.0);

        pocketfft::c2r(pocketfft::shape_t{out_data_real.size()},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
                       pocketfft::shape_t{0},
                       pocketfft::BACKWARD,
                       hilbert_time_buffer2.data(),
                       out_data_imag.data(),
                       1.0);
    }

    void update_with_new_phase() {
        if(use_averaging) {
            phase_count++;
            for(size_t i = 0; i < num_freqs; i++) {
                double delta = spectral_phase[i] - spectral_phase_mean[i];
                spectral_phase_mean[i] += delta / phase_count;
                double delta2 = spectral_phase[i] - spectral_phase_mean[i];
                spectral_phase_mean_sqr[i] += delta * delta2;
            }
        } else {
            phase_count = 0;
            for(size_t i = 0; i < num_freqs; i++) {
                spectral_phase_mean[i] = spectral_phase[i];
                spectral_phase_mean_sqr[i] = 0.0;
            }
        }
    }

    void calculate_group_delay() {
        double dw_factor = 0.5 / (d_f * twoPi<double>());
        group_delay[0] = 2 * dw_factor * (spectral_phase_mean[1] - spectral_phase_mean[0]);
        for(size_t i = 1; i < num_freqs - 1; i++) {
            group_delay[i] = dw_factor * (spectral_phase_mean[i + 1] - spectral_phase_mean[i - 1]);
        }
        group_delay[num_freqs - 1] =
            2 * dw_factor *
            (spectral_phase_mean[num_freqs - 1] - spectral_phase_mean[num_freqs - 2]);
    }

  public:
    SpectralInterferometry() { setup_fft(); }
    ~SpectralInterferometry() {}
    void set_averaging(bool new_state) { use_averaging = new_state; }
    void reset_frequencies(size_t n, double f_min, double f_max) {
        if(n < 2)
            return;
        if(num_freqs == n && min_f == f_min && max_f == f_max)
            return;
        frequencies = std::vector<double>(n);
        d_f = (f_max - f_min) / (n - 1);
        num_freqs = n;
        min_f = f_min;
        max_f = f_max;
        for(size_t i = 0; i < n; i++) {
            frequencies[i] = f_min + static_cast<double>(i) * d_f;
        }
        if(fft_size != num_freqs) {
            setup_fft();
        }
        if(reference_data_a.size() > 0) {
            reference_data_a_interpolated =
                wavelength_to_frequency(frequencies, wavelengths, reference_data_a);
        }
        if(reference_data_b.size() > 0) {
            reference_data_b_interpolated =
                wavelength_to_frequency(frequencies, wavelengths, reference_data_b);
        }

        interference_data = std::vector<double>(wavelengths);
        interference_data_interpolated = std::vector<double>(num_freqs, 0.0);
        spectral_phase = std::vector<double>(num_freqs, 0.0);
        spectral_phase_mean = spectral_phase;
        spectral_phase_mean_sqr = spectral_phase;
        group_delay = spectral_phase;
        phase_count = 0;
        is_configured = true;
    }

    void reset_phase() {
        phase_count = 0;
        spectral_phase_mean = std::vector<double>(num_freqs, 0.0);
        spectral_phase_mean_sqr = spectral_phase_mean;
    }

    void acquire_new_phase(Spectrometer &s) {
        interference_data_interpolated = s.acquire_single_frequency(frequencies);
        calculate_phase();
        update_with_new_phase();
        calculate_group_delay();
    }

    void acquire_new_interferogram(Spectrometer &s) {
        interference_data_interpolated = s.acquire_single_frequency(frequencies);
    }
    void calculate_phase() {
        for(size_t i = 0; i < num_freqs; i++) {
            hilbert_real_buffer[i] = interference_data_interpolated[i] -
                                     reference_data_b_interpolated[i] -
                                     reference_data_a_interpolated[i];
        }
        filtered_hilbert(hilbert_real_buffer,
                         hilbert_real_buffer,
                         hilbert_imag_buffer);
        spectral_phase = std::vector<double>(num_freqs);
        for(size_t i = 0; i < num_freqs; i++) {
            spectral_phase[i] = std::atan2(hilbert_real_buffer[i], hilbert_imag_buffer[i]);
        }
        unwrap(spectral_phase);
    }
    double *get_phase_data() { return spectral_phase.data(); }
    double *get_mean_phase_data() { return spectral_phase_mean.data(); }
    double *get_group_delay_data() { return group_delay.data(); }
    double *get_frequencies() { return frequencies.data(); }
    size_t get_num_freqs() { return num_freqs; }
    size_t get_num_times() { return num_freqs / 2 + 1; }
    double *get_interference_data() { return interference_data_interpolated.data(); }
    double *get_reference_a() { return reference_data_a_interpolated.data(); }
    double *get_reference_b() { return reference_data_b_interpolated.data(); }
    std::vector<double> &get_reference_a_vector() { return reference_data_a_interpolated; }
    std::vector<double> &get_reference_b_vector() { return reference_data_b_interpolated; }

    double *get_time_data() { return fft_data_real.data(); }

    double *get_time_reference_a() { return fft_reference_a_real.data(); }

    double *get_time_reference_b() { return fft_reference_b_real.data(); }
    double *get_time_scale() { return fft_scale.data(); }
    double *get_time_filter() { return time_filter.data(); }
    void set_time_filter(double t0, double sigma, double ord) {
        filter_t0 = t0;
        filter_width = sigma;
        filter_order = ord;
    }
    void generate_time_plot() {
        pocketfft::r2c(pocketfft::shape_t{interference_data_interpolated.size()},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
                       pocketfft::shape_t{0},
                       pocketfft::FORWARD,
                       interference_data_interpolated.data(),
                       hilbert_time_buffer.data(),
                       1.0);

        pocketfft::r2c(pocketfft::shape_t{reference_data_a_interpolated.size()},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
                       pocketfft::shape_t{0},
                       pocketfft::FORWARD,
                       reference_data_a_interpolated.data(),
                       fft_reference_a.data(),
                       1.0);

        pocketfft::r2c(pocketfft::shape_t{reference_data_b_interpolated.size()},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
                       pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
                       pocketfft::shape_t{0},
                       pocketfft::FORWARD,
                       reference_data_b_interpolated.data(),
                       fft_reference_b.data(),
                       1.0);
        double dt = 1.0e3 / (num_freqs * d_f);
        double max_signal = 0.0;
        for(size_t i = 0; i < (num_freqs / 2 + 1); i++) {
            constexpr auto MODULUS_SQUARED = [](const std::complex<double> &x) {
                return x.real() * x.real() + x.imag() * x.imag();
            };
            fft_reference_a_real[i] = MODULUS_SQUARED(fft_reference_a[i]);
            fft_reference_b_real[i] = MODULUS_SQUARED(fft_reference_b[i]);
            fft_data_real[i] =
                MODULUS_SQUARED(hilbert_time_buffer[i] - fft_reference_a[i] - fft_reference_b[i]);
            max_signal = maxN(fft_reference_a_real[i], max_signal);
            time_filter[i] = std::exp(
                -std::pow((static_cast<double>(i) * dt - filter_t0) / filter_width, filter_order) /
                sqrt(2.0));
            fft_scale[i] = static_cast<double>(i) * dt;
        }
        for(size_t i = 0; i < (num_freqs / 2 + 1); i++) {
            time_filter[i] *= max_signal;
        }
    }

    bool check_configuration_status() { return is_configured; }
    void acquire_reference_a(BatchAcquisition &batch_control,
                             const size_t n,
                             const double integration_time,
                             const double seconds_to_wait,
                             Spectrometer &s) {
        if(n == 0)
            return;
        batch_control.acquire_batch(n, integration_time, seconds_to_wait, s);
        wavelengths = s.wavelengths_copy();
        reference_data_a = std::vector<double>(wavelengths.size(), 0.0);
        double norm_factor = 1.0 / static_cast<double>(n);
        std::vector<double> full_set = batch_control.get_data_vector();
        for(size_t i = 0; i < wavelengths.size(); i++) {
            for(size_t j = 0; j < n; j++) {
                reference_data_a[i] += full_set[i + wavelengths.size() * j];
            }
            reference_data_a[i] *= norm_factor;
        }
        reference_data_a_interpolated =
            wavelength_to_frequency(frequencies, wavelengths, reference_data_a);
    }

    void acquire_reference_b(BatchAcquisition &batch_control,
                             const size_t n,
                             const double integration_time,
                             const double seconds_to_wait,
                             Spectrometer &s) {
        if(n == 0)
            return;
        batch_control.acquire_batch(n, integration_time, seconds_to_wait, s);
        wavelengths = s.wavelengths_copy();
        reference_data_b = std::vector<double>(wavelengths.size(), 0.0);
        double norm_factor = 1.0 / static_cast<double>(n);
        std::vector<double> full_set = batch_control.get_data_vector();
        for(size_t i = 0; i < wavelengths.size(); i++) {
            for(size_t j = 0; j < n; j++) {
                reference_data_b[i] += full_set[i + j * wavelengths.size()];
            }
            reference_data_b[i] *= norm_factor;
        }
        reference_data_b_interpolated =
            wavelength_to_frequency(frequencies, wavelengths, reference_data_b);
    }

    void save(std::string &path, bool timestamp) {
        if(!is_configured)
            return;
        if(timestamp) {
            auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
                                 std::chrono::system_clock::now().time_since_epoch())
                                 .count();
            size_t last_period = path.find_last_of(".");
            // if the extension is weird or absent, just put timestamp at end
            if(last_period == std::string::npos || (path.length() - last_period) > 6) {
                path.append(std::format("{}", timestamp));
            } else {
                path.insert(last_period, std::format("{}", timestamp));
            }
        }
        std::ofstream fs(path, std::ios::binary);
        if(fs.fail())
            return;
        fs.precision(10);
        calculate_group_delay();
        for(size_t j = 0; j < num_freqs; j++) {
            fs << frequencies[j];
            fs << " ";
            fs << spectral_phase_mean[j];
            fs << " ";
            fs << group_delay[j];
            fs << " ";
            fs << reference_data_a_interpolated[j];
            fs << " ";
            fs << reference_data_b_interpolated[j];
            fs << " ";
            fs << interference_data_interpolated[j];
            fs << '\x0A';
        }
    }
};

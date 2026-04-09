#pragma once
#include <vector>
#include <chrono>
#include <complex>
#include <memory>
#include <span>
#include <format>
#include <fstream>
#include "ExternalLibraries/pocketfft_hdronly.h"
#include "spectrometer.hpp"
#include "batch_acquisition.hpp"

class SpectralInterferometry {
    std::vector<double> wavelengths;
    std::vector<double> frequencies;
    std::vector<double> wavelengths_on_frequency_grid;
    std::vector<double> interference_data;
    std::vector<double> interference_data_interpolated;
    std::vector<double> reference_data_A;
    std::vector<double> reference_data_A_interpolated;
    std::vector<double> reference_data_B;
    std::vector<double> reference_data_B_interpolated;

    std::vector<double> spectral_phase;
    std::vector<double> spectral_phase_mean;
    std::vector<double> spectral_phase_mean_sqr;
    std::vector<double> group_delay;
    size_t phase_count = 0;


    std::vector<std::complex<double>> fft_data;
    std::vector<std::complex<double>> fft_reference_A;
    std::vector<std::complex<double>> fft_reference_B;
    std::vector<double> fft_data_real;
    std::vector<double> fft_reference_A_real;
    std::vector<double> fft_reference_B_real;
    std::vector<double> fft_scale;
    std::vector<double> time_filter;
    std::vector<std::complex<double>> hilbert_time_buffer;
    std::vector<std::complex<double>> hilbert_time_buffer2;
    std::vector<double> hilbert_real_buffer;
    std::vector<double> hilbert_imag_buffer;
    bool has_data = false;
    bool is_configured = false;
    bool use_averaging = true;
    double dF = 2e12;
    double minF = 100e12;
    double maxF = 2e12 * 512 + 100e12;
    double filterT0 = 90;
    double filter_width = 50;
    double filter_order = 8;
    size_t num_freqs = 512;
    size_t fft_size = 0;
    void unwrap(std::vector<double>& phaseInOut) {
        double offset{};
        for (int i = 1; i < num_freqs; i++) {
            double delta = phaseInOut[i] - phaseInOut[i - 1];
            delta -= twoPi<double>() * round(delta / twoPi<double>());
            phaseInOut[i] = phaseInOut[i - 1] + delta;
        }
    }
    void setup_fft() {
        std::unique_ptr<double[]> setupWorkD(new double[num_freqs]);
        std::unique_ptr<std::complex<double>[]> setupWorkC(new std::complex<double>[num_freqs]);
        fft_size = num_freqs;
        hilbert_time_buffer = std::vector<std::complex<double>>(num_freqs/2 + 1);
        hilbert_time_buffer2 = hilbert_time_buffer;
        hilbert_real_buffer = std::vector<double>(num_freqs);
        hilbert_imag_buffer = hilbert_real_buffer;
        fft_reference_A = hilbert_time_buffer;
        fft_reference_B = hilbert_time_buffer;
        time_filter = std::vector<double>(num_freqs / 2 + 1);
        fft_scale = time_filter;
        fft_reference_A_real = time_filter;
        fft_reference_B_real = time_filter;
        fft_data_real = time_filter;
    }


    void filtered_hilbert(std::vector<double>& inData,
        std::vector<double>& outDataReal,
        std::vector<double>& outDataImag,
        double filter0,
        double filterSigma,
        int filterOrd) {
        if (fft_size != num_freqs) {
            setup_fft();
        }

        pocketfft::r2c(
            pocketfft::shape_t{inData.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::shape_t{0},
            pocketfft::FORWARD,
            inData.data(),
            hilbert_time_buffer.data(),
            1.0);

        double dt = 1.0 / (num_freqs * dF);
        std::complex<double> ii(0.0, 1.0);
        for (int i = 0; i < (num_freqs / 2 + 1); i++) {
            hilbert_time_buffer[i] *=
                std::exp(-std::pow(
                    (static_cast<double>(i) * dt - filter0) / filterSigma, filterOrd)
                    / sqrt(2.0));
            hilbert_time_buffer2[i] = ii * hilbert_time_buffer[i];
        }

        pocketfft::c2r(
            pocketfft::shape_t{outDataReal.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::shape_t{0},
            pocketfft::BACKWARD,
            hilbert_time_buffer.data(),
            outDataReal.data(),
            1.0);

        pocketfft::c2r(
            pocketfft::shape_t{outDataReal.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::shape_t{0},
            pocketfft::BACKWARD,
            hilbert_time_buffer2.data(),
            outDataImag.data(),
            1.0);
    }

    void update_with_new_phase() {
        if (use_averaging) {
            phase_count++;
            for (int i = 0; i < num_freqs; i++) {
                double delta = spectral_phase[i] - spectral_phase_mean[i];
                spectral_phase_mean[i] += delta / phase_count;
                double delta2 = spectral_phase[i] - spectral_phase_mean[i];
                spectral_phase_mean_sqr[i] += delta * delta2;
            }
        }
        else {
            phase_count = 0;
            for (int i = 0; i < num_freqs; i++) {
                spectral_phase_mean[i] = spectral_phase[i];
                spectral_phase_mean_sqr[i] = 0.0;
            }
        }
    }

    void calculate_group_delay() {
        double dwFactor = 0.5/(dF * twoPi<double>());
        group_delay[0] = 2 * dwFactor * (spectral_phase_mean[1] - spectral_phase_mean[0]);
        for (int i = 1; i < num_freqs-1; i++) {
            group_delay[i] = dwFactor * (spectral_phase_mean[i + 1] - spectral_phase_mean[i - 1]);
        }
        group_delay[num_freqs - 1] = 2 * dwFactor * (spectral_phase_mean[num_freqs - 1] - spectral_phase_mean[num_freqs - 2]);
    }

public:
    SpectralInterferometry() {
        setup_fft();
    }
    ~SpectralInterferometry() {
    }
    void set_averaging(bool newState) {
        use_averaging = newState;
    }
    void reset_frequencies(size_t N, double fMin, double fMax) {
        if (N < 2) return;
        if (num_freqs == N && minF == fMin && maxF == fMax) return;
        frequencies = std::vector<double>(N);
        dF = (fMax - fMin) / (N - 1);
        num_freqs = N;
        minF = fMin;
        maxF = fMax;
        for (int i = 0; i < N; i++) {
            frequencies[i] = fMin + static_cast<double>(i) * dF;
        }
        if (fft_size != num_freqs) {
            setup_fft();
        }
        if (reference_data_A.size() > 0) {
            reference_data_A_interpolated = wavelength_to_frequency(frequencies, wavelengths, reference_data_A);
        }
        if (reference_data_B.size() > 0) {
            reference_data_B_interpolated = wavelength_to_frequency(frequencies, wavelengths, reference_data_B);
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

    void acquire_new_phase(Spectrometer& s) {
        interference_data_interpolated = s.acquire_single_frequency(frequencies);
        calculate_phase();
        update_with_new_phase();
        calculate_group_delay();
    }

    void acquire_new_interferogram(Spectrometer& s) {
        interference_data_interpolated = s.acquire_single_frequency(frequencies);
    }
    void calculate_phase() {
        for (int i = 0; i < num_freqs; i++) {
            hilbert_real_buffer[i] = interference_data_interpolated[i] - reference_data_B_interpolated[i] - reference_data_A_interpolated[i];
        }
        filtered_hilbert(hilbert_real_buffer, hilbert_real_buffer, hilbert_imag_buffer, filterT0, filter_width, 8);
        spectral_phase = std::vector<double>(num_freqs);
        for (int i = 0; i < num_freqs; i++) {
            spectral_phase[i] = std::atan2(hilbert_real_buffer[i], hilbert_imag_buffer[i]);
        }
        unwrap(spectral_phase);
    }
    double* get_phase_data() {
        return spectral_phase.data();
    }
    double* get_mean_phase_data() {
        return spectral_phase_mean.data();
    }
    double* get_group_delay_data() {
        return group_delay.data();
    }
    double* get_frequencies() {
        return frequencies.data();
    }
    size_t get_num_freqs() {
        return num_freqs;
    }
    size_t get_num_times() {
        return num_freqs / 2 + 1;
    }
    double* get_interference_data() {
        return interference_data_interpolated.data();
    }
    double* get_reference_A() {
        return reference_data_A_interpolated.data();
    }
    double* get_reference_B() {
        return reference_data_B_interpolated.data();
    }
    std::vector<double>& get_reference_A_vector() {
        return reference_data_A_interpolated;
    }
    std::vector<double>& get_reference_B_vector() {
        return reference_data_B_interpolated;
    }

    double* get_time_data() {
        return fft_data_real.data();
    }

    double* get_time_reference_A() {
        return fft_reference_A_real.data();
    }

    double* get_time_reference_B() {
        return fft_reference_B_real.data();
    }
    double* get_time_scale() {
        return fft_scale.data();
    }
    double* get_time_filter() {
        return time_filter.data();
    }
    void set_time_filter(double t0, double sigma, double ord) {
        filterT0 = t0;
        filter_width = sigma;
        filter_order = ord;
    }
    void generate_time_plot() {
        pocketfft::r2c(
            pocketfft::shape_t{interference_data_interpolated.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::shape_t{0},
            pocketfft::FORWARD,
            interference_data_interpolated.data(),
            hilbert_time_buffer.data(),
            1.0);

        pocketfft::r2c(
            pocketfft::shape_t{reference_data_A_interpolated.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::shape_t{0},
            pocketfft::FORWARD,
            reference_data_A_interpolated.data(),
            fft_reference_A.data(),
            1.0);

        pocketfft::r2c(
            pocketfft::shape_t{reference_data_B_interpolated.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::shape_t{0},
            pocketfft::FORWARD,
            reference_data_B_interpolated.data(),
            fft_reference_B.data(),
            1.0);
        double dt = 1.0e3 / (num_freqs * dF);
        double max_signal = 0.0;
        for (int i = 0; i < (num_freqs / 2 + 1); i++) {
            constexpr auto modulus_squared = [](const std::complex<double>& x){
                return x.real() * x.real() + x.imag() * x.imag();
            };
            fft_reference_A_real[i] = modulus_squared(fft_reference_A[i]);
            fft_reference_B_real[i] = modulus_squared(fft_reference_B[i]);
            fft_data_real[i] = modulus_squared(hilbert_time_buffer[i] - fft_reference_A[i] - fft_reference_B[i]);
            max_signal = maxN(fft_reference_A_real[i], max_signal);
            time_filter[i] = std::exp(-std::pow(
                (static_cast<double>(i) * dt - filterT0) / filter_width, filter_order)
                / sqrt(2.0));
            fft_scale[i] = static_cast<double>(i) * dt;
        }
        for (int i = 0; i < (num_freqs / 2 + 1); i++) {
            time_filter[i] *= max_signal;
        }
    }

    bool check_configuration_status() {
        return is_configured;
    }
    void acquire_reference_A(BatchAcquisition& batchControl, const size_t N, const double integration_time, const double seconds_to_wait, Spectrometer& s) {
        if (N == 0) return;
        batchControl.acquire_batch(N, integration_time, seconds_to_wait, s);
        wavelengths = s.wavelengths_copy();
        reference_data_A = std::vector<double>(wavelengths.size(), 0.0);
        double normFactor = 1.0 / static_cast<double>(N);
        std::vector<double> fullSet = batchControl.get_data_vector();
        for (int i = 0; i < wavelengths.size(); i++) {
            for (int j = 0; j < N; j++) {
                reference_data_A[i] += fullSet[i + wavelengths.size() * j];
            }
            reference_data_A[i] *= normFactor;
        }
        reference_data_A_interpolated = wavelength_to_frequency(frequencies, wavelengths, reference_data_A);
    }

    void acquireReferenceB(BatchAcquisition& batchControl, const size_t N, const double integration_time, const double seconds_to_wait, Spectrometer& s) {
        if (N == 0) return;
        batchControl.acquire_batch(N, integration_time, seconds_to_wait, s);
        wavelengths = s.wavelengths_copy();
        reference_data_B = std::vector<double>(wavelengths.size(), 0.0);
        double normFactor = 1.0 / static_cast<double>(N);
        std::vector<double> fullSet = batchControl.get_data_vector();
        for (int i = 0; i < wavelengths.size(); i++) {
            for (int j = 0; j < N; j++) {
                reference_data_B[i] +=fullSet[i + j *wavelengths.size()];
            }
            reference_data_B[i] *= normFactor;
        }
        reference_data_B_interpolated = wavelength_to_frequency(frequencies, wavelengths, reference_data_B);
    }

    void save(std::string& path, bool timestamp) {
        if (timestamp) {
            auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            size_t lastPeriod = path.find_last_of(".");
            //if the extension is weird or absent, just put timestamp at end
            if (lastPeriod == std::string::npos || (path.length() - lastPeriod) > 6) {
                path.append(std::format("{}", timestamp));
            }
            else {
                path.insert(lastPeriod, std::format("{}", timestamp));
            }
        }
        std::ofstream fs(path, std::ios::binary);
        if (fs.fail()) return;
        fs.precision(10);
        for (size_t j = 0; j < num_freqs; j++) {
            fs << frequencies[j];
            fs << " ";
            fs << spectral_phase_mean[j];
            fs << " ";
            fs << std::sqrt(spectral_phase_mean_sqr[j] / (phase_count - 1));
            fs << " ";
            fs << reference_data_A_interpolated[j];
            fs << " ";
            fs << reference_data_B_interpolated[j];
            fs << " ";
            fs << interference_data_interpolated[j];
            fs << '\x0A';
        }
    }
};

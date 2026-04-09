#include <iostream>
#include <thread>
#include <chrono>
#include <algorithm>
#include <complex>
#include <memory>
#include <span>
#include "ExternalLibraries/pocketfft_hdronly.h"
#include "api/OceanDirectAPI.h"
#include "ExternalLibraries/LightwaveExplorerGraphicalClasses.h"

void destroyMainWindowCallback();
bool updateDisplay();
void handleRunButton();
void drawSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectrumFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectraFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferenceSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferenceSpectrumTime(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferencePhase(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferencegroup_delay(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drop_down_change_callback();
void handleget_overlay0();
void handleget_overlay1();
void handleget_overlay2();
void handleCollapsePanel();
void handleDeleteOverlay0();
void handleDeleteOverlay1();
void handleDeleteOverlay2();
void handleGetdark_spectrum();
void handleDeletedark_spectrum();
void handleRefreshRequest();
void savePathCallback();
void handleSave();
void handleReferenceA();
void handleReferenceB();
void handleResetController();
void handlereset_phase();
void handleSavePhase();
void svgCallback();

double constexpr modulus_squared(const std::complex<double>& x){
    return x.real() * x.real() + x.imag() * x.imag();
}

//c++ port of the Rust code I wrote for the attoworld library
constexpr std::vector<double> fornberg_stencil(const size_t order, const std::span<const double>& positions, const double position_out){
    size_t n_pos = positions.size();
    size_t cols = order + 1;
    std::vector<double> delta_current(n_pos * cols);
    std::vector<double> delta_previous(n_pos * cols);
    delta_current[0] = 1.0;
    double c1 = 1.0;
    for (size_t n=1; n<n_pos; n++) {
        delta_previous.swap(delta_current);
        double c2 = 1.0;
        size_t zero_previous = n <= order;
        size_t min_n_order = std::min(n, order);
        for( size_t v = 0; v<n; v++) {
            double c3 = positions[n] - positions[v];
            c2 *= c3;

            if (zero_previous) {
                delta_previous[n * n_pos + v] = 0.0;
            }

            for (size_t m = 0; m<=min_n_order; m++) {
                double last_element = m == 0 ?
                    0.0 :
                    static_cast<double>(m) * delta_previous[(m - 1) * n_pos + v];

                delta_current[m * n_pos + v] = ((positions[n] - position_out)
                    * delta_previous[m * n_pos + v]
                    - last_element)
                    / c3;
            }
        }

        for (size_t m = 0; m<=min_n_order; m++) {
            double first_element = m == 0 ?
                0.0 :
                static_cast<double>(m) * delta_previous[(m - 1) * n_pos + n - 1];

            delta_current[m * n_pos + n] = (c1 / c2)
                * (first_element
                    - (positions[n - 1] - position_out) * delta_previous[m * n_pos + n - 1]);
        }

        c1 = c2;
    }
    return std::vector<double>(delta_current.begin() + order*n_pos, delta_current.begin() + cols*n_pos);
}

constexpr inline double interp_fornberg(const double target_frequency, const std::span<const double> frequencies, const std::span<const double> wavelengths, const std::span<const double> spectrum_in){
    const size_t stencil_width = 4;
    const size_t stencil_offset = stencil_width/2u;
    const double target_wavelength = constProd(lightC<double>(), 1e-3) / target_frequency;
    size_t stencil_x_start = std::distance(wavelengths.begin(), std::lower_bound(wavelengths.begin(), wavelengths.end(), target_wavelength));
    if(stencil_x_start < stencil_offset) {
        stencil_x_start = 0;
    }
    else if(stencil_x_start >= wavelengths.size() - stencil_width){
        stencil_x_start = wavelengths.size() - stencil_width;
    }
    else{
        if(abs(wavelengths[stencil_x_start] - target_wavelength) > abs(wavelengths[stencil_x_start + 1] - target_wavelength)){
            stencil_x_start++;
        }
        stencil_x_start -= stencil_offset;
    }
    std::span<const double> stencil_wavelengths = std::span(wavelengths).subspan(stencil_x_start, stencil_width);
    auto stencil = fornberg_stencil(0, stencil_wavelengths, target_wavelength);
    double interpolated_value = 0.0;
    for(size_t i = 0; i < stencil_width; i++){
        interpolated_value += stencil[i] * spectrum_in[stencil_x_start + i];
    }
    return target_wavelength * target_wavelength * 1e-6 * interpolated_value;
}

std::vector<double> wavelength_to_frequency(const std::span<const double> frequencies, const std::span<const double> wavelengths, const std::span<const double> spectrum_in)
{
    std::vector<double> spectrum_out(frequencies.size());
#pragma omp parallel for
    for (int j = 0; j < frequencies.size(); j++) {
        spectrum_out[j] = interp_fornberg(frequencies[j], frequencies, wavelengths, spectrum_in);
    }
    return spectrum_out;
}


//Spectrometer class which handles the acquisition from the spectrometer, as well as the associated
//buffers for the data and overlays
class Spectrometer {
public:
    long device_id;
    int error;
    int pixel_count;
    std::vector<double> read_buffer;
    std::vector<double> read_buffer_minus_dark;
    std::vector<double> wavelengths_buffer;
    std::vector<double> overlay0;
    std::vector<double> overlay0_minus_dark;
    std::vector<double> overlay1;
    std::vector<double> overlay1_minus_dark;
    std::vector<double> overlay2;
    std::vector<double> overlay2_minus_dark;
    std::vector<double> overlay0F;
    std::vector<double> overlay1F;
    std::vector<double> overlay2F;
    std::vector<double> dark_spectrum;
    bool hasdark_spectrum = false;
    int lastOverlay = -1;
    bool isInitialized = false;
    bool isLocked = false;
    Spectrometer(long _device_id, int _pixel_count) :
        device_id(_device_id),
        pixel_count(_pixel_count),
        read_buffer(pixel_count),
        read_buffer_minus_dark(pixel_count),
        overlay0(pixel_count),
        overlay0_minus_dark(pixel_count),
        overlay1(pixel_count),
        overlay1_minus_dark(pixel_count),
        overlay2(pixel_count),
        overlay2_minus_dark(pixel_count),
        dark_spectrum(pixel_count),
        wavelengths_buffer(pixel_count){}
    void subtract_dark(std::vector<double>& dataVector, std::vector<double>& dataMinusDark) {
        if (!hasdark_spectrum) {
            dataMinusDark = dataVector;
            return;
        }

        for (size_t i = 0; i < pixel_count; i++) {
            dataMinusDark[i] = dataVector[i] - dark_spectrum[i];
        }
    }

    bool hasOverlay0 = false;
    bool hasOverlay1 = false;
    bool hasOverlay2 = false;
    virtual ~Spectrometer(){}
    virtual void set_integration_time(unsigned long integrationTimeMicroseconds) = 0;
    virtual void acquire_single() = 0;
    virtual std::vector<double> acquire_single_frequency(const std::vector<double>& frequencies) = 0;
    virtual void acquire_overlay(int overlayIndex) = 0;
    virtual void acquire_dark_spectrum() = 0;

    bool initialized() {
        return isInitialized;
    }
    void lock() {
        isLocked = true;
    }
    void unlock() {
        isLocked = false;
    }
    bool checkLock() {
        return isLocked;
    }
    void disabledark_spectrum() {
        hasdark_spectrum = false;
    }
    int get_overlay_count() {
        int overlayCount = 0;
        if (hasOverlay0) overlayCount++;
        if (hasOverlay1) overlayCount++;
        if (hasOverlay2) overlayCount++;
        return overlayCount;
    }
    double* get_overlay(int overlayIndex) {
        switch (overlayIndex) {
        case 0:
            if(hasdark_spectrum) return overlay0_minus_dark.data();
            return overlay0.data();
        case 1:
            if (hasdark_spectrum) return overlay1_minus_dark.data();
            return overlay1.data();
        case 2:
            if (hasdark_spectrum) return overlay2_minus_dark.data();
            return overlay2.data();
        default:
            return nullptr;
        }
    }

    double* get_overlayFrequency(int overlayIndex, const std::vector<double> frequencies) {
        switch (overlayIndex) {
        case 0:
            if (hasdark_spectrum) overlay0F = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay0_minus_dark);
            else overlay0F = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay0);
            return overlay0F.data();
        case 1:
            if (hasdark_spectrum) overlay1F = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay1_minus_dark);
            overlay1F = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay1);
            return overlay1F.data();
        case 2:
            if (hasdark_spectrum) overlay2F = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay2_minus_dark);
            overlay2F = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay2);
            return overlay2F.data();
        default:
            return nullptr;
        }
    }

    void deleteOverlay(int overlayIndex) {
        switch (overlayIndex) {
        case 0:
            hasOverlay0 = false;
            break;
        case 1:
            hasOverlay1 = false;
            break;
        case 2:
            hasOverlay2 = false;
            break;
        }
    }

    double* data() {
        if (hasdark_spectrum) { return read_buffer_minus_dark.data(); }
        else {
            return read_buffer.data();
        }
    }

    void appendBufferTo(std::vector<double>& outputBuffer) {
        if (!hasdark_spectrum) {
            outputBuffer.insert(outputBuffer.end(), read_buffer.begin(), read_buffer.end());
        }
        else {
            outputBuffer.insert(outputBuffer.end(), read_buffer_minus_dark.begin(), read_buffer_minus_dark.end());
        }
    }

    double* wavelengths() {
        return wavelengths_buffer.data();
    }

    std::vector<double> wavelengthsCopy() {
        return wavelengths_buffer;
    }

    size_t size() {
        return pixel_count;
    }

    int getErrorCode() {
        return error;
    }
};


class OceanSpectrometer : public Spectrometer{
public:
	OceanSpectrometer(long device_idinput) :
	    Spectrometer(device_idinput, odapi_get_formatted_spectrum_length(device_idinput, &error))
	{
        isInitialized = (error==0);
        if(isInitialized) odapi_get_wavelengths(device_id, &error, wavelengths_buffer.data(), pixel_count);
	}
    virtual ~OceanSpectrometer() override {
        odapi_close_device(device_id, &error);
    }

    void set_integration_time(unsigned long integrationTimeMicroseconds) override {
        odapi_set_integration_time_micros(device_id, &error, integrationTimeMicroseconds);
    }

    void acquire_single() override {
        odapi_get_formatted_spectrum(device_id, &error, read_buffer.data(), pixel_count);
        subtract_dark(read_buffer, read_buffer_minus_dark);
    }

    std::vector<double> acquire_single_frequency(const std::vector<double>& frequencies) override {
        odapi_get_formatted_spectrum(device_id, &error, read_buffer.data(), pixel_count);
        subtract_dark(read_buffer, read_buffer_minus_dark);
        return wavelength_to_frequency(frequencies, wavelengths_buffer, read_buffer_minus_dark);
    }

    void acquire_overlay(int overlayIndex) override {
        switch (overlayIndex) {
        case 0:
            odapi_get_formatted_spectrum(device_id, &error, overlay0.data(), pixel_count);
            subtract_dark(overlay0, overlay0_minus_dark);
            lastOverlay = 0;
            hasOverlay0 = true;
            break;
        case 1:
            odapi_get_formatted_spectrum(device_id, &error, overlay1.data(), pixel_count);
            subtract_dark(overlay1, overlay1_minus_dark);
            lastOverlay = 1;
            hasOverlay1 = true;
            break;
        case 2:
            odapi_get_formatted_spectrum(device_id, &error, overlay2.data(), pixel_count);
            subtract_dark(overlay2, overlay2_minus_dark);
            lastOverlay = 2;
            hasOverlay2 = true;
            break;
        }
    }
    void acquire_dark_spectrum() override {
        dark_spectrum = std::vector<double>(pixel_count);
        odapi_get_formatted_spectrum(device_id, &error, dark_spectrum.data(), pixel_count);
        hasdark_spectrum = true;
    }
};

std::vector<std::unique_ptr<Spectrometer>> spectrometerSet;

//Batch acquisition class which holds the buffers for acquiring a series of shots, plus the methods
//to acquire and save the data
class BatchAcquisition {
    std::vector<double> data;
    std::vector<double> wavelengths;
    bool has_data = false;
    bool acquisition_finished = false;
    size_t spectrum_size = 0;
    size_t acquisition_count = 0;
public:
    void acquire_batch(const size_t N, const double integrationTime, const double secondsToWait, Spectrometer& s) {
        if (N == 0) return;
        s.lock();
        spectrum_size = s.size();
        wavelengths = s.wavelengthsCopy();
        data.clear();
        data.reserve(N * spectrum_size);
        acquisition_count = 0;
        acquisition_finished = false;
        s.set_integration_time((unsigned long)round(1000 * integrationTime));
        for (size_t i = 0; i < N; i++) {
            s.acquire_single();
            s.appendBufferTo(data);
            acquisition_count++;
            std::this_thread::sleep_for(std::chrono::milliseconds((size_t)(1000.0 * secondsToWait)));
            has_data = true;
        }
        acquisition_finished = true;
        s.unlock();
    }

    double* get_data(size_t offset) {
        return &data.data()[offset * spectrum_size];
    }

    std::vector<double>& get_data_vector() {
        return data;
    }

    void save(std::string& path, bool timestamp) {
        if (!has_data) return;
        if (timestamp) {
            auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            size_t lastPeriod = path.find_last_of(".");
            //if the extension is weird or absent, just put timestamp at end
            if (lastPeriod == std::string::npos || (path.length() - lastPeriod) > 6) {
                path.append(Sformat("{}", timestamp));
            }
            else {
                path.insert(lastPeriod, Sformat("{}", timestamp));
            }
        }
        std::ofstream fs(path, std::ios::binary);
        if (fs.fail()) return;
        fs.precision(10);
        for (size_t j = 0; j < spectrum_size; j++) {
            fs << wavelengths[j];
            for (size_t i = 0; i < acquisition_count; i++) {
                fs << " ";
                fs << data[i * spectrum_size + j];
            }
            fs << '\x0A';
        }
    }

    size_t get_spectrum_size() {
        return spectrum_size;
    }
};
BatchAcquisition theBatch;

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


    void filteredHilbert(std::vector<double>& inData,
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
    void setAveraging(bool newState) {
        use_averaging = newState;
    }
    void resetFrequencies(size_t N, double fMin, double fMax) {
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
        filteredHilbert(hilbert_real_buffer, hilbert_real_buffer, hilbert_imag_buffer, filterT0, filter_width, 8);
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
        //fftw_execute_dft_r2c(fftPlanD2Z, reference_data_A_interpolated.data(), (fftw_complex*)fft_reference_A.data());
        pocketfft::r2c(
            pocketfft::shape_t{reference_data_A_interpolated.size()},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(double))},
            pocketfft::stride_t{static_cast<ptrdiff_t>(sizeof(std::complex<double>))},
            pocketfft::shape_t{0},
            pocketfft::FORWARD,
            reference_data_A_interpolated.data(),
            fft_reference_A.data(),
            1.0);
        //fftw_execute_dft_r2c(fftPlanD2Z, reference_data_A_interpolated.data(), (fftw_complex*)fft_reference_B.data());
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
    void acquire_reference_A(BatchAcquisition& batchControl, const size_t N, const double integrationTime, const double secondsToWait, Spectrometer& s) {
        if (N == 0) return;
        batchControl.acquire_batch(N, integrationTime, secondsToWait, s);
        wavelengths = s.wavelengthsCopy();
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

    void acquireReferenceB(BatchAcquisition& batchControl, const size_t N, const double integrationTime, const double secondsToWait, Spectrometer& s) {
        if (N == 0) return;
        batchControl.acquire_batch(N, integrationTime, secondsToWait, s);
        wavelengths = s.wavelengthsCopy();
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
                path.append(Sformat("{}", timestamp));
            }
            else {
                path.insert(lastPeriod, Sformat("{}", timestamp));
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
SpectralInterferometry theInterferenceController;

//Main class for controlling the interface
class MainGui {
    bool queueUpdate = false;
    bool queueSliderUpdate = false;
    bool queueSliderMove = false;
    bool queueInterfaceValuesUpdate = false;
    bool queueSavePathUpdate = false;
    bool isRunningLive = false;
    bool noSpectrometers = false;
    int sliderTarget = 0;

public:
    LweTextBox textBoxes[54];
    LweButton buttons[18];
    LweButton miniButtons[12];
    LweConsole console;
    LweConsole sequence;
    LweConsole fitCommand;
    LweTextBox filePaths[4];
    LwePulldown pulldowns[10];
    LweDrawBox drawBoxes[8];
    LweDrawBox progressBarBox;
    LweCheckBox checkBoxes[4];
    LweSlider plotSlider;
    LweWindow window;
    std::string pathBuffer;
    int saveSVG = 0;
    bool loadedDefaults = false;
    unsigned int timeoutID = 0;
    void requestPlotUpdate() {
        queueUpdate = true;
    }

    void requestLive() {
        isRunningLive = true;
    }

    void stopLive() {
        isRunningLive = false;
    }
    bool noSpectrometersFound() {
        return noSpectrometers;
    }
    [[nodiscard]] constexpr bool runningLive() { return isRunningLive; }
    void applyUpdate() {
        if (queueUpdate) {
            queueUpdate = false;
            drawBoxes[0].queueDraw();
        }
        if (queueSavePathUpdate && pathBuffer != "?LWE_LOADING??") {
            queueSavePathUpdate = false;
            if (pathBuffer == "?LWE_NOPATH??") return;
            filePaths[0].overwritePrint(pathBuffer);
        }
    }
    void requestSliderUpdate() {
        queueSliderUpdate = true;
    }
    void requestSliderMove(int target) {
        queueSliderMove = true;
        sliderTarget = target;
    }
    void requestInterfaceValuesUpdate() {
        queueInterfaceValuesUpdate = true;
    }
    void requestSavePathUpdate() {
#ifndef __APPLE__
        pathBuffer = std::string("?LWE_LOADING??");
        pathFromSaveDialog(pathBuffer, "txt", "Text (.txt)");
#endif
        queueSavePathUpdate = true;
    }
    void activate(GtkApplication* app) {
        int buttonWidth = 4;
        int smallButton = 3;
        int textWidth = 2;
        int labelWidth = 2;
        int plotWidth = 12;
        int plotHeight = 6;
        int pathChars = 42;
        int colWidth = labelWidth + 2 * textWidth;
        int textCol0 = 0;
        int textCol1 = textWidth;
        int textCol2 = 2 * textWidth;
        int textCol3 = 3 * textWidth;
        int textCol4 = 4 * textWidth;
        int textCol5 = 5 * textWidth;
        int textCol6 = 6 * textWidth;
        int buttonCol1 = 8;

        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", true, NULL);
        window.init(app, "Spectrometer Control and Recording for Attosecond Beamlines", 1400, 800);
        GtkWidget* parentHandle = window.parentHandle();

        buttons[0].init(("Run"), parentHandle, buttonCol1, 1, buttonWidth, 1, handleRunButton);
        buttons[1].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 3, buttonWidth/2, 1, handleget_overlay0);
        buttons[2].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 4, buttonWidth / 2, 1, handleget_overlay1);
        buttons[3].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 5, buttonWidth / 2, 1, handleget_overlay2);
        buttons[3].init(("Acquire"), parentHandle, buttonCol1, 7, buttonWidth, 1, handleSave);

        buttons[4].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 3, buttonWidth / 2, 1, handleDeleteOverlay0);
        buttons[5].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 4, buttonWidth / 2, 1, handleDeleteOverlay1);
        buttons[6].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 5, buttonWidth / 2, 1, handleDeleteOverlay2);
        buttons[7].init(("\xf0\x9f\x95\xaf\xef\xb8\x8f"), parentHandle, buttonCol1, 2, buttonWidth / 2, 1, handleGetdark_spectrum);
        buttons[8].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1 + buttonWidth / 2, 2, buttonWidth / 2, 1, handleDeletedark_spectrum);

        buttons[9].init(("Ref. A"), parentHandle, 0, 12, smallButton, 1, handleReferenceA);
        buttons[10].init(("Ref. B"), parentHandle, smallButton, 12, smallButton, 1, handleReferenceB);
        buttons[11].init(("Reset"), parentHandle, smallButton * 2, 12, smallButton, 1, handlereset_phase);
        buttons[12].init(("Save"), parentHandle, smallButton * 3, 12, smallButton, 1, handleSavePhase);

        //RGB active
        textBoxes[1].init(parentHandle, textCol1, 2, textWidth, 1);
        textBoxes[2].init(parentHandle, textCol2, 2, textWidth, 1);
        textBoxes[3].init(parentHandle, textCol3, 2, textWidth, 1);

        //RGB overlay0
        textBoxes[4].init(parentHandle, textCol1, 3, textWidth, 1);
        textBoxes[5].init(parentHandle, textCol2, 3, textWidth, 1);
        textBoxes[6].init(parentHandle, textCol3, 3, textWidth, 1);

        //RGB overlay1
        textBoxes[7].init(parentHandle, textCol1, 4, textWidth, 1);
        textBoxes[8].init(parentHandle, textCol2, 4, textWidth, 1);
        textBoxes[9].init(parentHandle, textCol3, 4, textWidth, 1);

        //RGB overlay2
        textBoxes[10].init(parentHandle, textCol1, 5, textWidth, 1);
        textBoxes[11].init(parentHandle, textCol2, 5, textWidth, 1);
        textBoxes[12].init(parentHandle, textCol3, 5, textWidth, 1);

        textBoxes[0].init(parentHandle, textCol3, 7, textWidth, 1);
        textBoxes[0].setLabel(-3 * textWidth, 0, "Exposure (ms)");
        textBoxes[0].overwritePrint(std::string("40"));

        textBoxes[13].init(parentHandle, textCol3, 8, textWidth, 1);
        textBoxes[13].setLabel(-3 * textWidth, 0, "Pause, count (s, #)");
        textBoxes[14].init(parentHandle, textCol4, 8, textWidth, 1);
        textBoxes[13].overwritePrint(std::string("0"));
        textBoxes[14].overwritePrint(std::string("10"));


        textBoxes[15].init(parentHandle, textCol3, 9, 2, 1);
        textBoxes[15].setMaxCharacters(6);
        textBoxes[16].init(parentHandle, textCol4, 9, 2, 1);
        textBoxes[16].setMaxCharacters(6);
        textBoxes[17].init(parentHandle, textCol5, 9, 2, 1);
        textBoxes[17].setMaxCharacters(6);
        textBoxes[15].setLabel(-3 * textWidth, 0, "Freqs. (#, min,max)");
        textBoxes[15].overwritePrint(std::string("2048"));

        textBoxes[18].init(parentHandle, textCol3, 11, 2, 1);
        textBoxes[18].setMaxCharacters(6);
        textBoxes[19].init(parentHandle, textCol4, 11, 2, 1);
        textBoxes[19].setMaxCharacters(6);
        textBoxes[20].init(parentHandle, textCol5, 11, 2, 1);
        textBoxes[20].setMaxCharacters(6);
        textBoxes[18].setLabel(-3 * textWidth, 0, "Filter (t0, sig., ord.)");
        textBoxes[18].overwritePrint(std::string("300"));
        textBoxes[19].overwritePrint(std::string("80"));
        textBoxes[20].overwritePrint(std::string("12"));

        filePaths[0].init(parentHandle, 0, 6, 10, 1);
        filePaths[0].setMaxCharacters(pathChars);
        filePaths[0].overwritePrint(Sformat("DefaultOutput.txt"));
        checkBoxes[2].init("\xe2\x8c\x9a", parentHandle, buttonCol1 + buttonWidth/2, 8, 2, 1);
        buttons[7].init(("..."), parentHandle, buttonCol1 + buttonWidth / 2, 6, buttonWidth / 2, 1, savePathCallback, 0);
        textBoxes[48].init(window.parentHandle(4), 3, 0, 2, 1);
        textBoxes[49].init(window.parentHandle(4), 5, 0, 2, 1);
        textBoxes[50].init(window.parentHandle(4), 9, 0, 2, 1);
        textBoxes[51].init(window.parentHandle(4), 11, 0, 2, 1);
        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawSpectrum);
        checkBoxes[1].init(("Log"), window.parentHandle(4), 13, 0, 1, 1);
        checkBoxes[3].init(("Average phase"), window.parentHandle(4), 15, 0, 1, 1);
        buttons[13].init(("xlim"), window.parentHandle(4), 2, 0, 1, 1, handleRefreshRequest);
        buttons[13].setTooltip("Apply the entered x limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[13].squeeze();
        buttons[14].init(("ylim"), window.parentHandle(4), 7, 0, 1, 1, handleRefreshRequest);
        buttons[14].setTooltip("Apply the entered y limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[14].squeeze();
        buttons[15].init(("SVG"), window.parentHandle(4), 1, 0, 1, 1, svgCallback);
        buttons[15].setTooltip("Generate SVG files of the four line plots, with filenames based on the base path set above");
        buttons[15].squeeze();
        buttons[16].init("\xe2\x86\x94\xef\xb8\x8f", window.parentHandle(4), 0, 0, 1, 1, handleCollapsePanel);
        buttons[16].setTooltip("Collapse/expand the data entry panel");

        for (int i = 0; i < 18; i++) {
            textBoxes[i].setMaxCharacters(6);
        }
        pulldowns[1].addElement("Spectrometer live (nm)");
        pulldowns[1].addElement("Spectrometer live (THz)");
        pulldowns[1].addElement("All spectrometers");
        pulldowns[1].addElement("Interferometry: Spectrum");
        pulldowns[1].addElement("Interferometry: Time");
        pulldowns[1].addElement("Interferometry: phase");
        pulldowns[1].addElement("Interferometry: Group delay");
        pulldowns[1].init(parentHandle, 0, 1, 8, 1);
        g_signal_connect(pulldowns[1].elementHandle, "notify::selected", G_CALLBACK(drop_down_change_callback), NULL);
        gtk_widget_set_visible(buttons[9].elementHandle, false);
        gtk_widget_set_visible(buttons[10].elementHandle, false);
        gtk_widget_set_visible(buttons[11].elementHandle, false);
        gtk_widget_set_visible(buttons[12].elementHandle, false);
        gtk_widget_set_visible(textBoxes[18].label, false);
        gtk_widget_set_visible(textBoxes[18].elementHandle, false);
        gtk_widget_set_visible(textBoxes[19].elementHandle, false);
        gtk_widget_set_visible(textBoxes[20].elementHandle, false);
        console.init(window.parentHandle(1), 0, 0, 1, 1);
        console.cPrint("Attached spectrometers:\n");

        window.present();
        initializeSpectrometers();
        pulldowns[0].init(parentHandle, 0, 0, 12, 1);
        drawBoxes[0].queueDraw();
        timeoutID = g_timeout_add(20, G_SOURCE_FUNC(updateDisplay), NULL);
    }

    void initializeSpectrometers() {
        int deviceCount = odapi_probe_devices();
        if (deviceCount == 0) {
            console.cPrint("I didn't find any spectrometers...\nMake sure you've installed the driver.\n");
            noSpectrometers = true;
            return;
        }
        int device_idCount = odapi_get_number_of_device_ids();

        std::vector<long> device_ids(device_idCount);
        int retrievedIdCount = odapi_get_device_ids(device_ids.data(), device_idCount);
        int error = 0;
        console.cPrint("RID count {}\n", retrievedIdCount);

        spectrometerSet = std::vector<std::unique_ptr<Spectrometer>>();

		for (int i = 0; i < device_idCount; i++) {
			odapi_open_device(device_ids[i], &error);
            if (error) {
                console.cPrint("Couldn't open the device! Error code{}\n", error);
                noSpectrometers = true;
                return;
            }
			// Get the device name
			const int nameLength = 128;
			char deviceName[nameLength] = { 0 };
			odapi_get_device_name(device_ids[i], &error, deviceName, nameLength);
			if (error != 0) {
				console.cPrint("Failed to retrieve the spectrometer type. The error code is:  {}\n", error);
                noSpectrometers = true;
                return;
			}
			// and serial number
			int serialNumberLength = odapi_get_serial_number_maximum_length(device_ids[i], &error);
            if (error != 0) {
                console.cPrint("Failed to retrieve the serial number length. The error code is:  {}\n", error);
                noSpectrometers = true;
                return;
            }
            std::unique_ptr<char> serialNumber(new char[serialNumberLength + 1] {});
			odapi_get_serial_number(device_ids[i], &error, serialNumber.get(), serialNumberLength);
			if (error != 0) {
				console.cPrint("Failed to retrieve the spectrometer serial number. The error code is:  {}\n", error);
                noSpectrometers = true;
                return;
			}
			else {
				console.cPrint("Device {}: {}\n    serial number: {}\n", i, deviceName, serialNumber.get());
                //console.cPrint("Device {}: {}\n    serial number: {}\n", i, deviceName, 0);
                std::string newElement = Sformat("{}: {}", i, deviceName);
                pulldowns[0].addElement(newElement.c_str());
			}
            spectrometerSet.push_back(std::make_unique<OceanSpectrometer>(device_ids[i]));
		}
	}

};
MainGui theGui;

void destroyMainWindowCallback() {
}

void handleRunButton() {
    if (theGui.runningLive()) {
        theGui.stopLive();
    }
    else {
        theGui.requestLive();
    }
}

void svgCallback() {
    theGui.saveSVG = 1;
    theGui.requestPlotUpdate();
}

void handleRefreshRequest() {
    theGui.requestPlotUpdate();
}

void handleget_overlay0() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_overlay(0);
}

void drop_down_change_callback(){
    bool visibility = theGui.pulldowns[1].getValue()>2;
    gtk_widget_set_visible(theGui.buttons[9].elementHandle, visibility);
    gtk_widget_set_visible(theGui.buttons[10].elementHandle, visibility);
    gtk_widget_set_visible(theGui.buttons[11].elementHandle, visibility);
    gtk_widget_set_visible(theGui.buttons[12].elementHandle, visibility);
    gtk_widget_set_visible(theGui.textBoxes[18].label, visibility);
    gtk_widget_set_visible(theGui.textBoxes[18].elementHandle, visibility);
    gtk_widget_set_visible(theGui.textBoxes[19].elementHandle, visibility);
    gtk_widget_set_visible(theGui.textBoxes[20].elementHandle, visibility);
}

void handleget_overlay1() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_overlay(1);
}

void handleget_overlay2() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_overlay(2);
}

void handleCollapsePanel() {
    theGui.window.toggleSettingsPanel();
}

void handleDeleteOverlay0() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).deleteOverlay(0);
}

void handleDeleteOverlay1() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).deleteOverlay(1);
}

void handleDeleteOverlay2() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).deleteOverlay(2);
}

void handleGetdark_spectrum() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_dark_spectrum();
}

void handleDeletedark_spectrum() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).disabledark_spectrum();
}

void acquisitionThread(int activeSpectrometer, size_t N, double integrationTime, double secondsToWait, std::string path, bool timestamp) {
    theBatch.acquire_batch(N, integrationTime, secondsToWait, (*spectrometerSet[activeSpectrometer]));
    theBatch.save(path, timestamp);
    theGui.console.tPrint("Finished writing {}!\n", path);
}

void handleSave() {
    theGui.stopLive();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((*spectrometerSet[activeSpectrometer]).checkLock()) return;
    std::string path;
    theGui.filePaths[0].copyBuffer(path);
    bool timestamp = theGui.checkBoxes[2].isChecked();
    size_t N = (size_t)theGui.textBoxes[14].valueDouble();
    double integrationTime = theGui.textBoxes[0].valueDouble();
    double waitTime = theGui.textBoxes[13].valueDouble();
    std::thread(acquisitionThread, activeSpectrometer, N, integrationTime, waitTime, path, timestamp).detach();
}

void handleSavePhase() {
    std::string path;
    theGui.filePaths[0].copyBuffer(path);
    bool timestamp = theGui.checkBoxes[2].isChecked();
    theInterferenceController.save(path, timestamp);
}

void referenceAAcquisitionThread(int activeSpectrometer, size_t N, double integrationTime, double secondsToWait) {
    theInterferenceController.acquire_reference_A(theBatch, N, integrationTime, secondsToWait, (*spectrometerSet[activeSpectrometer]));
}

void referenceBAcquisitionThread(int activeSpectrometer, size_t N, double integrationTime, double secondsToWait) {
    theInterferenceController.acquireReferenceB(theBatch, N, integrationTime, secondsToWait, (*spectrometerSet[activeSpectrometer]));
}

void handleReferenceA() {
    handleResetController();
    if (!theInterferenceController.check_configuration_status()) {
        theGui.console.cPrint("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stopLive();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((*spectrometerSet[activeSpectrometer]).checkLock()) return;
    size_t N = (size_t)theGui.textBoxes[14].valueDouble();
    double integrationTime = theGui.textBoxes[0].valueDouble();
    double waitTime = theGui.textBoxes[13].valueDouble();
    std::thread(referenceAAcquisitionThread, activeSpectrometer, N, integrationTime, waitTime).detach();
}

void handleReferenceB() {
    handleResetController();
    if (!theInterferenceController.check_configuration_status()) {
        theGui.console.cPrint("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stopLive();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((*spectrometerSet[activeSpectrometer]).checkLock()) return;
    size_t N = (size_t)theGui.textBoxes[14].valueDouble();
    double integrationTime = theGui.textBoxes[0].valueDouble();
    double waitTime = theGui.textBoxes[13].valueDouble();
    std::thread(referenceBAcquisitionThread, activeSpectrometer, N, integrationTime, waitTime).detach();
}

void handleResetController() {
    size_t num_freqs = static_cast<size_t>(theGui.textBoxes[15].valueDouble());
    if (num_freqs < 2) return;
    double fMin = theGui.textBoxes[16].valueDouble();
    double fMax = theGui.textBoxes[17].valueDouble();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (fMax == 0.0) {
        fMax = (*spectrometerSet[activeSpectrometer]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if (fMin == 0.0) {
        fMin = (*spectrometerSet[activeSpectrometer]).wavelengths()[(*spectrometerSet[activeSpectrometer]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    theInterferenceController.resetFrequencies(num_freqs, fMin, fMax);
}

void handlereset_phase() {
    theInterferenceController.reset_phase();
}

void drawSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    switch (theGui.pulldowns[1].getValue()) {
    case 0:
        break;
    case 1:
        drawSpectrumFrequency(area, cr, width, height, data);
        return;
    case 2:
        drawSpectraFrequency(area, cr, width, height, data);
        return;
    case 3:
        drawInterferenceSpectrum(area, cr, width, height, data);
        return;
    case 4:
        drawInterferenceSpectrumTime(area, cr, width, height, data);
        return;
    case 5:
        drawInterferencePhase(area, cr, width, height, data);
        return;
    case 6:
        drawInterferencegroup_delay(area, cr, width, height, data);
        return;
    default:
        return;
    }
    if (theGui.noSpectrometersFound()) return;
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).initialized()) {
        theGui.console.cPrint("Not initialized - error {}\n",(*spectrometerSet[0]).getErrorCode());
        return;
    }

    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;

    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceYmin = true;
        forceYmax = yMax != 0.0;
    }
    if(logPlot && !forceYmin){
        forceYmin = true;
        yMin = 1e-1;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_Spectrum.svg");
        sPlot.SVGPath = svgPath;
    }

    if (theGui.runningLive() && (activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).checkLock()) {
        (*spectrometerSet[activeSpectrometer]).set_integration_time((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        (*spectrometerSet[activeSpectrometer]).acquire_single();
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = (*spectrometerSet[activeSpectrometer]).wavelengths();
    sPlot.hasDataX = true;
    sPlot.data = (*spectrometerSet[activeSpectrometer]).data();
    sPlot.Npts = (*spectrometerSet[activeSpectrometer]).size();
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceYmin;
    sPlot.forceYmax = forceYmax;
    if (forceYmax)sPlot.forcedYmax = yMax;
    if (forceYmin)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Wavelength (nm)";
    sPlot.yLabel = "Spectrum (counts)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    if ((*spectrometerSet[activeSpectrometer]).get_overlay_count() > 0) {
        int firstAdded = -1;
        int secondAdded = -1;
        sPlot.ExtraLines = (*spectrometerSet[activeSpectrometer]).get_overlay_count();
        if ((*spectrometerSet[activeSpectrometer]).hasOverlay0) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay(0);
            firstAdded = 0;
        }
        else if ((*spectrometerSet[activeSpectrometer]).hasOverlay1) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay(1);
            firstAdded = 1;
        }
        else if ((*spectrometerSet[activeSpectrometer]).hasOverlay2) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay(2);
            firstAdded = 2;
        }
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && (*spectrometerSet[activeSpectrometer]).hasOverlay1) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlay(1);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && (*spectrometerSet[activeSpectrometer]).hasOverlay2) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlay(2);

        if (sPlot.ExtraLines > 2) sPlot.data4 = (*spectrometerSet[activeSpectrometer]).get_overlay(2);
    }
    sPlot.plot(cr);
}

void drawSpectrumFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {

    if (theGui.noSpectrometersFound()) return;
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).initialized()) {
        theGui.console.cPrint("Not initialized - error {}\n", (*spectrometerSet[0]).getErrorCode());
        return;
    }

    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;

    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceYmin = true;
        forceYmax = yMax != 0.0;
    }
    if(logPlot && !forceYmin){
        forceYmin = true;
        yMin = 1e-1;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t num_freqs = static_cast<size_t>(theGui.textBoxes[15].valueDouble());
    if (num_freqs < 2) return;
    double fMin = theGui.textBoxes[16].valueDouble();
    double fMax = theGui.textBoxes[17].valueDouble();

    if (fMax == 0.0) {
        fMax = (*spectrometerSet[activeSpectrometer]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if (fMin == 0.0) {
        fMin = (*spectrometerSet[activeSpectrometer]).wavelengths()[(*spectrometerSet[activeSpectrometer]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    double dF = (fMax - fMin) / static_cast<double>(num_freqs - 1);
    std::vector<double> frequencies(num_freqs);
    for (size_t i = 0; i < num_freqs; i++) {
        frequencies[i] = fMin + static_cast<double>(i) * dF;
    }
    std::vector<double> liveSpectrum;
    if (theGui.runningLive() && (activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).checkLock() && num_freqs > 0) {
        (*spectrometerSet[activeSpectrometer]).set_integration_time((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        liveSpectrum = (*spectrometerSet[activeSpectrometer]).acquire_single_frequency(frequencies);
    }
    else return;

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = frequencies.data();
    sPlot.hasDataX = true;
    sPlot.data = liveSpectrum.data();
    sPlot.Npts = num_freqs;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceYmin;
    sPlot.forceYmax = forceYmax;
    if (forceYmax)sPlot.forcedYmax = yMax;
    if (forceYmin)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Spectrum (counts)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    if ((*spectrometerSet[activeSpectrometer]).get_overlay_count() > 0) {
        int firstAdded = -1;
        int secondAdded = -1;
        sPlot.ExtraLines = (*spectrometerSet[activeSpectrometer]).get_overlay_count();
        if ((*spectrometerSet[activeSpectrometer]).hasOverlay0) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlayFrequency(0, frequencies);
            firstAdded = 0;
        }
        else if ((*spectrometerSet[activeSpectrometer]).hasOverlay1) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlayFrequency(1, frequencies);
            firstAdded = 1;
        }
        else if ((*spectrometerSet[activeSpectrometer]).hasOverlay2) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlayFrequency(2, frequencies);
            firstAdded = 2;
        }
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && (*spectrometerSet[activeSpectrometer]).hasOverlay1) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlayFrequency(1, frequencies);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && (*spectrometerSet[activeSpectrometer]).hasOverlay2) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlayFrequency(2, frequencies);

        if (sPlot.ExtraLines > 2) sPlot.data4 = (*spectrometerSet[activeSpectrometer]).get_overlayFrequency(2, frequencies);
    }
    sPlot.plot(cr);
}

void drawSpectraFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (theGui.noSpectrometersFound()) return;
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;

    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceYmin = true;
        forceYmax = yMax != 0.0;
    }
    if(logPlot && !forceYmin){
        forceYmin = true;
        yMin = 1e-1;
    }
    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t num_freqs = static_cast<size_t>(theGui.textBoxes[15].valueDouble());
    if (num_freqs < 2) return;
    double fMin = theGui.textBoxes[16].valueDouble();
    double fMax = theGui.textBoxes[17].valueDouble();
    if (fMax == 0.0) {
        fMax = (*spectrometerSet[0]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
        for (int i = 1; i < spectrometerSet.size(); i++) {
            fMax = maxN(fMax, 1e-3 * lightC<double>() / (*spectrometerSet[i]).wavelengths()[0]);
        }
    }
    if (fMin == 0.0) {
        fMin = (*spectrometerSet[0]).wavelengths()[(*spectrometerSet[0]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
        for (int i = 1; i < spectrometerSet.size(); i++) {
            fMin = minN(fMin, 1e-3 * lightC<double>() / (*spectrometerSet[i]).wavelengths()[(*spectrometerSet[0]).size() - 1]);
        }
    }

    double dF = (fMax - fMin) / static_cast<double>(num_freqs - 1);
    std::vector<double> frequencies(num_freqs);
    for (size_t i = 0; i < num_freqs; i++) {
        frequencies[i] = fMin + static_cast<double>(i) * dF;
    }
    std::vector<std::vector<double>> liveSpectra(spectrometerSet.size());
    if (theGui.runningLive() && num_freqs > 0) {
        for (int i = 0; i < spectrometerSet.size(); i++) {
            liveSpectra[i] = (*spectrometerSet[i]).acquire_single_frequency(frequencies);
        }
    }
    else return;

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = frequencies.data();
    sPlot.hasDataX = true;
    sPlot.data = liveSpectra[0].data();
    sPlot.Npts = num_freqs;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceYmin;
    sPlot.forceYmax = forceYmax;
    if (forceYmax)sPlot.forcedYmax = yMax;
    if (forceYmin)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Spectrum (counts)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    sPlot.ExtraLines = static_cast<int>(liveSpectra.size())-1;
    if (liveSpectra.size() > 1) sPlot.data2 = liveSpectra[1].data();
    if (liveSpectra.size() > 2) sPlot.data3 = liveSpectra[2].data();
    if (liveSpectra.size() > 3) sPlot.data4 = liveSpectra[3].data();
    sPlot.plot(cr);
}

void drawInterferenceSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (theGui.noSpectrometersFound()) return;
    handleResetController();
    if (!theInterferenceController.check_configuration_status()) handleResetController();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_InterferenceF.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }

    if (theGui.runningLive() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.setAveraging(theGui.checkBoxes[3].isChecked());
        theInterferenceController.acquire_new_interferogram((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.get_frequencies();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.get_interference_data();
    sPlot.Npts = theInterferenceController.get_num_freqs();
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Spectrum (counts)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    sPlot.ExtraLines = 0;
    if (theInterferenceController.get_reference_A_vector().size() == theInterferenceController.get_num_freqs()) {
        sPlot.data2 = theInterferenceController.get_reference_A();
        sPlot.ExtraLines = 1;
        if (theInterferenceController.get_reference_B_vector().size() == theInterferenceController.get_num_freqs()) {
            sPlot.data3 = theInterferenceController.get_reference_B();
            sPlot.ExtraLines = 2;
        }
    }
    sPlot.plot(cr);
}

void drawInterferenceSpectrumTime(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (theGui.noSpectrometersFound()) return;
    handleResetController();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_InterferenceF.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }
    double t0 = theGui.textBoxes[18].valueDouble();
    double sigma = theGui.textBoxes[19].valueDouble();
    double ord = theGui.textBoxes[20].valueDouble();
    if (theGui.runningLive() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.setAveraging(theGui.checkBoxes[3].isChecked());
        theInterferenceController.acquire_new_interferogram((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }
    theInterferenceController.set_time_filter(t0, sigma, ord);
    theInterferenceController.generate_time_plot();

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.get_time_scale();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.get_time_data();
    sPlot.Npts = theInterferenceController.get_num_times();
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Spectrum (counts)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    sPlot.data2 = theInterferenceController.get_time_reference_A();
    sPlot.data3 = theInterferenceController.get_time_reference_B();
    sPlot.data4 = theInterferenceController.get_time_filter();
    sPlot.ExtraLines = 3;
    sPlot.plot(cr);
}


void drawInterferencePhase(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (theGui.noSpectrometersFound()) return;
    handleResetController();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_InterferencePhase.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }

    if (theGui.runningLive() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.setAveraging(theGui.checkBoxes[3].isChecked());
        theInterferenceController.acquire_new_phase((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.get_frequencies();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.get_mean_phase_data();
    sPlot.Npts = theInterferenceController.get_num_freqs();
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Phase (rad)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    sPlot.ExtraLines = 0;
    sPlot.plot(cr);
}

void drawInterferencegroup_delay(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (theGui.noSpectrometersFound()) return;
    handleResetController();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_InterferencePhase.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[1].valueDouble() + theGui.textBoxes[2].valueDouble() + theGui.textBoxes[3].valueDouble())) {
        mainColor = LweColor(theGui.textBoxes[1].valueDouble(), theGui.textBoxes[2].valueDouble(), theGui.textBoxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.textBoxes[4].valueDouble() + theGui.textBoxes[5].valueDouble() + theGui.textBoxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.textBoxes[4].valueDouble(), theGui.textBoxes[5].valueDouble(), theGui.textBoxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.textBoxes[7].valueDouble() + theGui.textBoxes[8].valueDouble() + theGui.textBoxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.textBoxes[7].valueDouble(), theGui.textBoxes[8].valueDouble(), theGui.textBoxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.textBoxes[10].valueDouble() + theGui.textBoxes[11].valueDouble() + theGui.textBoxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.textBoxes[10].valueDouble(), theGui.textBoxes[11].valueDouble(), theGui.textBoxes[12].valueDouble(), 1);
    }

    if (theGui.runningLive() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.setAveraging(theGui.checkBoxes[3].isChecked());
        theInterferenceController.acquire_new_phase((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.get_frequencies();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.get_group_delay_data();
    sPlot.Npts = theInterferenceController.get_num_freqs();
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = mainColor;
    sPlot.color2 = overLay0Color;
    sPlot.color3 = overLay1Color;
    sPlot.color4 = overLay2Color;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Group delay (ps)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    sPlot.ExtraLines = 0;
    sPlot.plot(cr);
}

bool updateDisplay() {
    theGui.requestPlotUpdate();
    theGui.console.updateFromBuffer();
    theGui.applyUpdate();
    return true;
}

void savePathCallback() {
#ifdef __APPLE__
    theGui.pathBuffer = pathFromAppleSaveDialog();
#endif
    theGui.requestSavePathUpdate();
}

static void activate(GtkApplication* app, gpointer user_data) {
#if defined __linux__ || defined __APPLE__
    setlocale(LC_NUMERIC, "en_US.UTF-8");
#else
    setlocale(LC_NUMERIC, "en_US");
#endif
    theGui.activate(app);
}

int main(int argc, char** argv) {
    GtkApplication* app = gtk_application_new("io.github.NickKarpowicz.Scarab", (GApplicationFlags)0);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    int status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    return status;
}

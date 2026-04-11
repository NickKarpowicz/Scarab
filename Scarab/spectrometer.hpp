#pragma once
#include "interpolation.hpp"
#include <string>
#include <vector>
// Spectrometer class which handles the acquisition from the spectrometer, as well as the associated
// buffers for the data and overlays
class Spectrometer {
  public:
    int pixel_count;
    std::string name;
    std::string serial_number;
    std::vector<double> read_buffer;
    std::vector<double> read_buffer_minus_dark;
    std::vector<double> wavelengths_buffer;
    std::vector<double> overlay0;
    std::vector<double> overlay0_minus_dark;
    std::vector<double> overlay1;
    std::vector<double> overlay1_minus_dark;
    std::vector<double> overlay2;
    std::vector<double> overlay2_minus_dark;
    std::vector<double> overlay0_f;
    std::vector<double> overlay1_f;
    std::vector<double> overlay2_f;
    std::vector<double> dark_spectrum;
    bool has_dark_spectrum = false;
    int last_overlay = -1;
    bool is_initialized = false;
    bool is_locked = false;
    Spectrometer(int pixel_count)
        : pixel_count(pixel_count), read_buffer(pixel_count), read_buffer_minus_dark(pixel_count),
          wavelengths_buffer(pixel_count), overlay0(pixel_count), overlay0_minus_dark(pixel_count),
          overlay1(pixel_count), overlay1_minus_dark(pixel_count), overlay2(pixel_count),
          overlay2_minus_dark(pixel_count), dark_spectrum(pixel_count) {}
    void subtract_dark(std::vector<double> &data_vector, std::vector<double> &data_minus_dark) {
        if(!has_dark_spectrum) {
            data_minus_dark = data_vector;
            return;
        }

        for(int i = 0; i < pixel_count; i++) {
            data_minus_dark[i] = data_vector[i] - dark_spectrum[i];
        }
    }

    bool has_overlay_0 = false;
    bool has_overlay_1 = false;
    bool has_overlay_2 = false;
    virtual ~Spectrometer() {}
    virtual void set_integration_time(unsigned long integration_time_microseconds) = 0;
    virtual void acquire_single() = 0;
    virtual std::vector<double>
    acquire_single_frequency(const std::vector<double> &frequencies) = 0;
    virtual void acquire_overlay(int overlay_index) = 0;
    virtual void acquire_dark_spectrum() = 0;

    bool initialized() { return is_initialized; }
    void lock() { is_locked = true; }
    void unlock() { is_locked = false; }
    bool check_lock() { return is_locked; }
    void disable_dark_spectrum() { has_dark_spectrum = false; }
    int get_overlay_count() {
        int overlay_count = 0;
        if(has_overlay_0)
            overlay_count++;
        if(has_overlay_1)
            overlay_count++;
        if(has_overlay_2)
            overlay_count++;
        return overlay_count;
    }
    double *get_overlay(int overlay_index) {
        switch(overlay_index) {
            case 0:
                if(has_dark_spectrum)
                    return overlay0_minus_dark.data();
                return overlay0.data();
            case 1:
                if(has_dark_spectrum)
                    return overlay1_minus_dark.data();
                return overlay1.data();
            case 2:
                if(has_dark_spectrum)
                    return overlay2_minus_dark.data();
                return overlay2.data();
            default:
                return nullptr;
        }
    }

    double *get_overlay_frequency(int overlay_index, const std::vector<double> frequencies) {
        switch(overlay_index) {
            case 0:
                if(has_dark_spectrum)
                    overlay0_f = wavelength_to_frequency(frequencies,
                                                         wavelengths_buffer,
                                                         overlay0_minus_dark);
                else
                    overlay0_f = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay0);
                return overlay0_f.data();
            case 1:
                if(has_dark_spectrum)
                    overlay1_f = wavelength_to_frequency(frequencies,
                                                         wavelengths_buffer,
                                                         overlay1_minus_dark);
                overlay1_f = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay1);
                return overlay1_f.data();
            case 2:
                if(has_dark_spectrum)
                    overlay2_f = wavelength_to_frequency(frequencies,
                                                         wavelengths_buffer,
                                                         overlay2_minus_dark);
                overlay2_f = wavelength_to_frequency(frequencies, wavelengths_buffer, overlay2);
                return overlay2_f.data();
            default:
                return nullptr;
        }
    }

    void delete_overlay(int overlay_index) {
        switch(overlay_index) {
            case 0:
                has_overlay_0 = false;
                break;
            case 1:
                has_overlay_1 = false;
                break;
            case 2:
                has_overlay_2 = false;
                break;
        }
    }

    double *data() {
        if(has_dark_spectrum) {
            return read_buffer_minus_dark.data();
        } else {
            return read_buffer.data();
        }
    }

    void append_buffer_to(std::vector<double> &output_buffer) {
        if(!has_dark_spectrum) {
            output_buffer.insert(output_buffer.end(), read_buffer.begin(), read_buffer.end());
        } else {
            output_buffer.insert(output_buffer.end(),
                                 read_buffer_minus_dark.begin(),
                                 read_buffer_minus_dark.end());
        }
    }

    double *wavelengths() { return wavelengths_buffer.data(); }

    std::vector<double> wavelengths_copy() { return wavelengths_buffer; }

    size_t size() { return pixel_count; }
};

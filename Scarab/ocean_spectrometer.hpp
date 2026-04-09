#pragma once
#include "spectrometer.hpp"
#include "api/OceanDirectAPI.h"
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

    void set_integration_time(unsigned long integration_timeMicroseconds) override {
        odapi_set_integration_time_micros(device_id, &error, integration_timeMicroseconds);
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
            has_overlay_0 = true;
            break;
        case 1:
            odapi_get_formatted_spectrum(device_id, &error, overlay1.data(), pixel_count);
            subtract_dark(overlay1, overlay1_minus_dark);
            lastOverlay = 1;
            has_overlay_1 = true;
            break;
        case 2:
            odapi_get_formatted_spectrum(device_id, &error, overlay2.data(), pixel_count);
            subtract_dark(overlay2, overlay2_minus_dark);
            lastOverlay = 2;
            has_overlay_2 = true;
            break;
        }
    }
    void acquire_dark_spectrum() override {
        dark_spectrum = std::vector<double>(pixel_count);
        odapi_get_formatted_spectrum(device_id, &error, dark_spectrum.data(), pixel_count);
        hasdark_spectrum = true;
    }
};

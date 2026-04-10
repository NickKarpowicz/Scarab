#pragma once
#include <format>
#include <string>
#include "spectrometer.hpp"
#include "api/OceanDirectAPI.h"

//wrapper for the open_device function to call when initializing the base Spectrometer class
[[nodiscard]] inline long open_ocean_spectrometer_wrapper(long device_id, int* error){
    odapi_open_device(device_id, error);
    return device_id;
}

class OceanSpectrometer : public Spectrometer{
public:
    long device_id;
    int error=0;
	OceanSpectrometer(long device_id_input) :
	    Spectrometer(odapi_get_formatted_spectrum_length(open_ocean_spectrometer_wrapper(device_id_input, &error), &error)),
		device_id(device_id_input)
	{
        is_initialized = (error==0);
        if(is_initialized) odapi_get_wavelengths(device_id, &error, wavelengths_buffer.data(), pixel_count);
        // Get the device name
        const int name_length = 128;
        char device_name[name_length] = { 0 };
        odapi_get_device_name(device_id, &error, device_name, name_length);
        name = std::string(device_name);
        // and serial number
        int serial_number_length = odapi_get_serial_number_maximum_length(device_id, &error);
        serial_number = std::string(serial_number_length,'\0');
        odapi_get_serial_number(device_id, &error, serial_number.data(), serial_number_length);
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
            last_overlay = 0;
            has_overlay_0 = true;
            break;
        case 1:
            odapi_get_formatted_spectrum(device_id, &error, overlay1.data(), pixel_count);
            subtract_dark(overlay1, overlay1_minus_dark);
            last_overlay = 1;
            has_overlay_1 = true;
            break;
        case 2:
            odapi_get_formatted_spectrum(device_id, &error, overlay2.data(), pixel_count);
            subtract_dark(overlay2, overlay2_minus_dark);
            last_overlay = 2;
            has_overlay_2 = true;
            break;
        }
    }
    void acquire_dark_spectrum() override {
        dark_spectrum = std::vector<double>(pixel_count);
        odapi_get_formatted_spectrum(device_id, &error, dark_spectrum.data(), pixel_count);
        has_dark_spectrum = true;
    }
};

inline std::string open_ocean_spectrometers(std::vector<std::unique_ptr<Spectrometer>>& spectrometer_set){
    std::string output;
    int device_count = odapi_probe_devices();
    if (device_count == 0) {
        output+="I didn't find any spectrometers...\nMake sure you've installed the driver.\n";
        return output;
    }
    int device_id_count = odapi_get_number_of_device_ids();

    std::vector<long> device_ids(device_id_count);
    int retrieved_id_count = odapi_get_device_ids(device_ids.data(), device_id_count);
    output += std::format("Ocean Optics spectrometer count {}\n", retrieved_id_count);
    for (int i = 0; i < device_id_count; i++) {
            spectrometer_set.push_back(std::make_unique<OceanSpectrometer>(device_ids[i]));
            if(((OceanSpectrometer*)spectrometer_set.back().get())->error){
                output+=std::format("Opening spectrometer {} failed with error {}\n",i,((OceanSpectrometer*)spectrometer_set.back().get())->error);
            }
    }
    return output;
};

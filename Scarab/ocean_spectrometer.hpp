#pragma once
#include <format>
#include <string>
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

inline std::string open_ocean_spectrometers(std::vector<std::unique_ptr<Spectrometer>>& spectrometer_set){
    std::string output;
    int deviceCount = odapi_probe_devices();
    if (deviceCount == 0) {
        output+="I didn't find any spectrometers...\nMake sure you've installed the driver.\n";
        return output;
    }
    int device_idCount = odapi_get_number_of_device_ids();

    std::vector<long> device_ids(device_idCount);
    int retrievedIdCount = odapi_get_device_ids(device_ids.data(), device_idCount);
    int error = 0;
    output += std::format("RID count {}\n", retrievedIdCount);
    for (int i = 0; i < device_idCount; i++) {
    odapi_open_device(device_ids[i], &error);
                if (error) {
                    output+=std::format("Couldn't open the device! Error code{}\n", error);
                    return output;
                }
    // Get the device name
    const int nameLength = 128;
    char deviceName[nameLength] = { 0 };
    odapi_get_device_name(device_ids[i], &error, deviceName, nameLength);
    if (error != 0) {
    output+=std::format("Failed to retrieve the spectrometer type. The error code is:  {}\n", error);
                    return output;
    }
    // and serial number
    int serialNumberLength = odapi_get_serial_number_maximum_length(device_ids[i], &error);
                if (error != 0) {
                    output+=std::format("Failed to retrieve the serial number length. The error code is:  {}\n", error);
                    return output;
                }
                std::unique_ptr<char> serialNumber(new char[serialNumberLength + 1] {});
    odapi_get_serial_number(device_ids[i], &error, serialNumber.get(), serialNumberLength);
    if (error != 0) {
    output += std::format("Failed to retrieve the spectrometer serial number. The error code is:  {}\n", error);
                    return output;
    }
            spectrometer_set.push_back(std::make_unique<OceanSpectrometer>(device_ids[i]));
            spectrometer_set.back()->name=std::string(deviceName);
            spectrometer_set.back()->serial_number=std::string(serialNumber.get());
    }
    return output;
};

#pragma once
#include "spectrometer.hpp"
#include <chrono>
#include <cmath>
#include <format>
#include <fstream>
#include <string>
#include <thread>
#include <vector>
// Batch acquisition class which holds the buffers for acquiring a series of shots, plus the methods
// to acquire and save the data
class BatchAcquisition {
    std::vector<double> data;
    std::vector<double> wavelengths;
    bool has_data = false;
    bool acquisition_finished = false;
    size_t spectrum_size = 0;
    size_t acquisition_count = 0;

  public:
    void acquire_batch(const size_t N,
                       const double integration_time,
                       const double seconds_to_wait,
                       Spectrometer &s) {
        if(N == 0)
            return;
        s.lock();
        spectrum_size = s.size();
        wavelengths = s.wavelengths_copy();
        data.clear();
        data.reserve(N * spectrum_size);
        acquisition_count = 0;
        acquisition_finished = false;
        s.set_integration_time((unsigned long)round(1000 * integration_time));
        for(size_t i = 0; i < N; i++) {
            s.acquire_single();
            s.append_buffer_to(data);
            acquisition_count++;
            std::this_thread::sleep_for(
                std::chrono::milliseconds((size_t)(1000.0 * seconds_to_wait)));
            has_data = true;
        }
        acquisition_finished = true;
        s.unlock();
    }

    double *get_data(size_t offset) { return &data.data()[offset * spectrum_size]; }

    std::vector<double> &get_data_vector() { return data; }

    void save(std::string &path, bool timestamp) {
        if(!has_data)
            return;
        if(timestamp) {
            auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
                                 std::chrono::system_clock::now().time_since_epoch())
                                 .count();
            size_t lastPeriod = path.find_last_of(".");
            // if the extension is weird or absent, just put timestamp at end
            if(lastPeriod == std::string::npos || (path.length() - lastPeriod) > 6) {
                path.append(std::format("{}", timestamp));
            } else {
                path.insert(lastPeriod, std::format("{}", timestamp));
            }
        }
        std::ofstream fs(path, std::ios::binary);
        if(fs.fail())
            return;
        fs.precision(10);
        for(size_t j = 0; j < spectrum_size; j++) {
            fs << wavelengths[j];
            for(size_t i = 0; i < acquisition_count; i++) {
                fs << " ";
                fs << data[i * spectrum_size + j];
            }
            fs << '\x0A';
        }
    }

    size_t get_spectrum_size() { return spectrum_size; }
};

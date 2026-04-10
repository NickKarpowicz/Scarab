#pragma once
#include <vector>
#include <span>
#include <cmath>
#include "ExternalLibraries/LightwaveExplorerHelpers.h"
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
        if(std::abs(wavelengths[stencil_x_start] - target_wavelength) > std::abs(wavelengths[stencil_x_start + 1] - target_wavelength)){
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

inline std::vector<double> wavelength_to_frequency(const std::span<const double> frequencies, const std::span<const double> wavelengths, const std::span<const double> spectrum_in)
{
    std::vector<double> spectrum_out(frequencies.size());
    int len = static_cast<int>(frequencies.size());
#pragma omp parallel for
    for (int j = 0; j < len; j++) {
        spectrum_out[j] = interp_fornberg(frequencies[j], frequencies, wavelengths, spectrum_in);
    }
    return spectrum_out;
}

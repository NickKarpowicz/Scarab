#include <iostream>
#include <thread>
#include <chrono>
#include <algorithm>
#include <complex>
#include <fftw3_mkl.h>
#include "api/OceanDirectAPI.h"
#include "LightwaveExplorerGTK/LightwaveExplorerGraphicalClasses.h"

void destroyMainWindowCallback();
bool updateDisplay();
void handleRunButton();
void drawSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectrumFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectraFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferenceSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferenceSpectrumTime(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferencePhase(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawInterferenceGroupDelay(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void handleGetOverlay0();
void handleGetOverlay1();
void handleGetOverlay2();
void handleDeleteOverlay0();
void handleDeleteOverlay1();
void handleDeleteOverlay2();
void handleGetDarkSpectrum();
void handleDeleteDarkSpectrum();
void handleRefreshRequest();
void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget);
void handleSave();
void handleReferenceA();
void handleReferenceB();
void handleResetController();
void handleResetPhase();
void handleSavePhase();
void svgCallback();

double modSquared(const std::complex<double>& x){
    return x.real() * x.real() + x.imag() * x.imag();
}

std::vector<double> wavelengthToFrequency(std::vector<double> frequencies, const std::vector<double>&wavelengths, const std::vector<double>&spectrumIn)
{
    std::vector<double> spectrumOut(frequencies.size());
    //lambda for performing single point interpolation
    auto interpolateSingle = [&](double targetFrequency) {
        double targetWavelength = constProd(lightC<double>(),1e-3) / targetFrequency;
        if (targetFrequency == 0.0 || targetWavelength < wavelengths[1] || targetWavelength > wavelengths[wavelengths.size() - 1]) {
            return 0.0;
        }
        // Find the index i such that wavelengths[i] <= targetWavelength <= wavelengths[i+1]
        size_t i = std::distance(wavelengths.begin(), std::lower_bound(wavelengths.begin(), wavelengths.end(), targetWavelength)) - 1;
        return (1e-6*targetWavelength*targetWavelength)*(spectrumIn[i - 1] + (spectrumIn[i] - spectrumIn[i - 1]) / (wavelengths[i] - wavelengths[i - 1]) * (targetWavelength - wavelengths[i - 1]));
    };

    for (size_t j = 0; j < frequencies.size(); j++) {
        spectrumOut[j] = interpolateSingle(frequencies[j]);
    }
    return spectrumOut;
}


//Spectrometer class which handles the acquisition from the spectrometer, as well as the associated
//buffers for the data and overlays
class OceanSpectrometer {
    std::vector<double> readBuffer;
    std::vector<double> readBufferMinusDark;
    std::vector<double> wavelengthsBuffer;
    std::vector<double> overlay0;
    std::vector<double> overlay0MinusDark;
    std::vector<double> overlay1;
    std::vector<double> overlay1MinusDark;
    std::vector<double> overlay2;
    std::vector<double> overlay2MinusDark;
    std::vector<double> overlay0F;
    std::vector<double> overlay1F;
    std::vector<double> overlay2F;
    std::vector<double> darkSpectrum;
    bool hasDarkSpectrum = false;
    int lastOverlay = -1;
    int pixelCount = 0;
    bool isInitialized = false;
    long deviceID = 0;
    int error = 0;
    bool isLocked = false;
    void subtractDark(std::vector<double>& dataVector, std::vector<double>& dataMinusDark) {
        if (!hasDarkSpectrum) {
            dataMinusDark = dataVector;
            return;
        }
        for (size_t i = 0; i < pixelCount; i++) {
            dataMinusDark[i] = dataVector[i] - darkSpectrum[i];
        }
    }
public:
    bool hasOverlay0 = false;
    bool hasOverlay1 = false;
    bool hasOverlay2 = false;
	void init(long deviceIDinput) {
		deviceID = deviceIDinput;
		pixelCount = odapi_get_formatted_spectrum_length(deviceID, &error);
        readBuffer = std::vector<double>(pixelCount);
        readBufferMinusDark = std::vector<double>(pixelCount);
        overlay0 = std::vector<double>(pixelCount);
        overlay0MinusDark = std::vector<double>(pixelCount);
        overlay1 = std::vector<double>(pixelCount);
        overlay1MinusDark = std::vector<double>(pixelCount);
        overlay2 = std::vector<double>(pixelCount);
        overlay2MinusDark = std::vector<double>(pixelCount);
        darkSpectrum = std::vector<double>(pixelCount);
        wavelengthsBuffer = std::vector<double>(pixelCount);
		
        if (error != 0 ) {
            isInitialized = false;
		}
		else {
            isInitialized = true;
            odapi_get_wavelengths(deviceID, &error, wavelengthsBuffer.data(), pixelCount);
		}
	}
    bool initialized() {
        return isInitialized;
    }

    void release() {
        odapi_close_device(deviceID, &error);
    }

    void setIntegrationTime(unsigned long integrationTimeMicroseconds) {
        odapi_set_integration_time_micros(deviceID, &error, integrationTimeMicroseconds);
    }

    void acquireSingle() {
        odapi_get_formatted_spectrum(deviceID, &error, readBuffer.data(), pixelCount); 
        subtractDark(readBuffer, readBufferMinusDark);
    }

    std::vector<double> acquireSingleFrequency(const std::vector<double>& frequencies) {
        odapi_get_formatted_spectrum(deviceID, &error, readBuffer.data(), pixelCount);
        subtractDark(readBuffer, readBufferMinusDark);
        return wavelengthToFrequency(frequencies, wavelengthsBuffer, readBuffer);
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
    void acquireOverlay(int overlayIndex) {        
        switch (overlayIndex) {
        case 0:
            odapi_get_formatted_spectrum(deviceID, &error, overlay0.data(), pixelCount);
            subtractDark(overlay0, overlay0MinusDark);
            lastOverlay = 0;
            hasOverlay0 = true;
            break;
        case 1:
            odapi_get_formatted_spectrum(deviceID, &error, overlay1.data(), pixelCount);
            subtractDark(overlay1, overlay1MinusDark);
            lastOverlay = 1;
            hasOverlay1 = true;
            break;
        case 2:
            odapi_get_formatted_spectrum(deviceID, &error, overlay2.data(), pixelCount);
            subtractDark(overlay2, overlay2MinusDark);
            lastOverlay = 2;
            hasOverlay2 = true;
            break;
        }
    }
    void acquireDarkSpectrum() {
        darkSpectrum = std::vector<double>(pixelCount);
        odapi_get_formatted_spectrum(deviceID, &error, darkSpectrum.data(), pixelCount);
        hasDarkSpectrum = true;
    }
    void disableDarkSpectrum() {
        hasDarkSpectrum = false;
    }

    int getOverlayCount() {
        int overlayCount = 0;
        if (hasOverlay0) overlayCount++;
        if (hasOverlay1) overlayCount++;
        if (hasOverlay2) overlayCount++;
        return overlayCount;
    }

    double* getOverlay(int overlayIndex) {
        switch (overlayIndex) {
        case 0:
            if(hasDarkSpectrum) return overlay0MinusDark.data();
            return overlay0.data();
        case 1:
            if (hasDarkSpectrum) return overlay1MinusDark.data();
            return overlay1.data();
        case 2:
            if (hasDarkSpectrum) return overlay2MinusDark.data();
            return overlay2.data();
        default:
            return nullptr;
        }
    }

    double* getOverlayFrequency(int overlayIndex, const std::vector<double> frequencies) {
        switch (overlayIndex) {
        case 0:
            if (hasDarkSpectrum) overlay0F = wavelengthToFrequency(frequencies, wavelengthsBuffer, overlay0MinusDark);
            else overlay0F = wavelengthToFrequency(frequencies, wavelengthsBuffer, overlay0);
            return overlay0F.data();
        case 1:
            if (hasDarkSpectrum) overlay1F = wavelengthToFrequency(frequencies, wavelengthsBuffer, overlay1MinusDark);
            overlay1F = wavelengthToFrequency(frequencies, wavelengthsBuffer, overlay1);
            return overlay1F.data();
        case 2:
            if (hasDarkSpectrum) overlay2F = wavelengthToFrequency(frequencies, wavelengthsBuffer, overlay2MinusDark);
            overlay2F = wavelengthToFrequency(frequencies, wavelengthsBuffer, overlay2);
            return overlay2F.data();
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
        if (hasDarkSpectrum) { return readBufferMinusDark.data(); }
        else {
            return readBuffer.data();
        }
        
    }

    void appendBufferTo(std::vector<double>& outputBuffer) {
        if (!hasDarkSpectrum) { 
            outputBuffer.insert(outputBuffer.end(), readBuffer.begin(), readBuffer.end()); 
        }
        else {
            outputBuffer.insert(outputBuffer.end(), readBufferMinusDark.begin(), readBufferMinusDark.end());
        }
        
    }

    double* wavelengths() {
        return wavelengthsBuffer.data();
    }

    std::vector<double> wavelengthsCopy() {
        return wavelengthsBuffer;
    }

    size_t size() {
        return pixelCount;
    }

    int getErrorCode() {
        return error;
    }
};
std::vector<OceanSpectrometer> spectrometerSet;

//Batch acquisition class which holds the buffers for acquiring a series of shots, plus the methods
//to acquire and save the data
class batchAcquisition {
    std::vector<double> data;
    std::vector<double> wavelengths;
    bool hasData = false;
    bool acquisitionFinished = false;
    size_t spectrumSize = 0;
    size_t acquisitionCount = 0;
public:

    void acquireBatch(const size_t N, const double integrationTime, const double secondsToWait, OceanSpectrometer& s) {
        if (N == 0) return;
        s.lock();
        spectrumSize = s.size();
        wavelengths = s.wavelengthsCopy();
        data.clear();
        data.reserve(N * spectrumSize);
        acquisitionCount = 0;
        acquisitionFinished = false;
        s.setIntegrationTime((unsigned long)round(1000 * integrationTime));
        for (size_t i = 0; i < N; i++) {
            s.acquireSingle();
            s.appendBufferTo(data);
            acquisitionCount++;
            std::this_thread::sleep_for(std::chrono::milliseconds((size_t)(1000.0 * secondsToWait)));
            hasData = true;
        }
        acquisitionFinished = true;
        s.unlock();
    }

    double* getData(size_t offset) {
        return &data.data()[offset * spectrumSize];
    }

    std::vector<double>& getDataVector() {
        return data;
    }

    void save(std::string& path, bool timestamp) {
        if (!hasData) return;
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
        for (size_t j = 0; j < spectrumSize; j++) {
            fs << wavelengths[j];
            for (size_t i = 0; i < acquisitionCount; i++) {
                fs << " ";
                fs << data[i * spectrumSize + j];
            }
            fs << '\x0A';
        }
    }

    size_t getSpectrumSize() {
        return spectrumSize;
    }
};
batchAcquisition theBatch;

class spectralInterferometry {
    std::vector<double> wavelengths;
    std::vector<double> frequencies;
    std::vector<double> wavelengthsOnFrequencyGrid;
    std::vector<double> interferenceData;
    std::vector<double> interferenceDataInterpolated;
    std::vector<double> referenceDataA;
    std::vector<double> referenceDataAInterpolated;
    std::vector<double> referenceDataB;
    std::vector<double> referenceDataBInterpolated;

    std::vector<double> spectralPhase;
    std::vector<double> spectralPhaseMean;
    std::vector<double> spectralPhaseM2;
    std::vector<double> groupDelay;
    size_t phaseCount = 0;


    std::vector<std::complex<double>> fftData;
    std::vector<std::complex<double>> fftReferenceA;
    std::vector<std::complex<double>> fftReferenceB;
    std::vector<double> fftDataReal;
    std::vector<double> fftReferenceAReal;
    std::vector<double> fftReferenceBReal;
    std::vector<double> fftScale;
    std::vector<double> timeFilter;
    std::vector<std::complex<double>> hilbertTimeBuffer;
    std::vector<std::complex<double>> hilbertTimeBuffer2;
    std::vector<double> hilbertRealBuffer;
    std::vector<double> hilbertImagBuffer;
    fftw_plan fftPlanD2Z;
    fftw_plan fftPlanZ2D;
    bool hasData = false;
    bool isConfigured = false;
    double dF = 2e12;
    double minF = 100e12;
    double maxF = 2e12 * 512 + 100e12;
    double filterT0 = 90;
    double filterWidth = 50;
    double filterOrder = 8;
    size_t Nfreq = 512;
    size_t FFTsize = 0;
    void unwrap(std::vector<double>& phaseInOut) {
        double offset{};
        for (int i = 1; i < Nfreq; i++) {
            double delta = phaseInOut[i] - phaseInOut[i - 1];
            delta -= twoPi<double>() * round(delta / twoPi<double>());
            phaseInOut[i] = phaseInOut[i - 1] + delta;
        }
    }
    void setupFFT() {
        std::unique_ptr<double[]> setupWorkD(new double[Nfreq]);
        std::unique_ptr<std::complex<double>[]> setupWorkC(new std::complex<double>[Nfreq]);
        const int fftw1[1] = { static_cast<int>(Nfreq) };
        fftPlanD2Z = fftw_plan_many_dft_r2c(1, fftw1, 1, setupWorkD.get(), fftw1, 1, static_cast<int>(Nfreq), (fftw_complex*)setupWorkC.get(), fftw1, 1, static_cast<int>(Nfreq / 2 + 1), FFTW_MEASURE);
        fftPlanZ2D = fftw_plan_many_dft_c2r(1, fftw1, 1, (fftw_complex*)setupWorkC.get(), fftw1, 1, static_cast<int>(Nfreq / 2 + 1), setupWorkD.get(), fftw1, 1, static_cast<int>(Nfreq), FFTW_MEASURE);
        FFTsize = Nfreq;
        hilbertTimeBuffer = std::vector<std::complex<double>>(Nfreq/2 + 1);
        hilbertTimeBuffer2 = hilbertTimeBuffer;
        hilbertRealBuffer = std::vector<double>(Nfreq);
        hilbertImagBuffer = hilbertRealBuffer;
        fftReferenceA = hilbertTimeBuffer;
        fftReferenceB = hilbertTimeBuffer;
        timeFilter = std::vector<double>(Nfreq / 2 + 1);
        fftScale = timeFilter;
        fftReferenceAReal = timeFilter;
        fftReferenceBReal = timeFilter;
        fftDataReal = timeFilter;
        
    }
    void destroyFFT() {
        fftw_destroy_plan(fftPlanD2Z);
        fftw_destroy_plan(fftPlanZ2D);
    }

    void filteredHilbert(std::vector<double>& inData,
        std::vector<double>& outDataReal,
        std::vector<double>& outDataImag,
        double filter0,
        double filterSigma,
        int filterOrd) {
        if (FFTsize != Nfreq) {
            destroyFFT();
            setupFFT();
        }
        fftw_execute_dft_r2c(fftPlanD2Z, inData.data(), (fftw_complex*)hilbertTimeBuffer.data());
        double dt = 1.0 / (Nfreq * dF);
        std::complex<double> ii(0.0, 1.0);
        for (int i = 0; i < (Nfreq / 2 + 1); i++) {
            hilbertTimeBuffer[i] *= 
                std::exp(-std::pow(
                    (static_cast<double>(i) * dt - filter0) / filterSigma, filterOrd) 
                    / sqrt(2.0));
            hilbertTimeBuffer2[i] = ii * hilbertTimeBuffer[i];
        }
        fftw_execute_dft_c2r(fftPlanZ2D, (fftw_complex*)hilbertTimeBuffer.data(), outDataReal.data());
        fftw_execute_dft_c2r(fftPlanZ2D, (fftw_complex*)hilbertTimeBuffer2.data(), outDataImag.data());
    }

    void updateWithNewPhase() {
        phaseCount++;
        for (int i = 0; i < Nfreq; i++) {
            double delta = spectralPhase[i] - spectralPhaseMean[i];
            spectralPhaseMean[i] += delta / phaseCount;
            double delta2 = spectralPhase[i] - spectralPhaseMean[i];
            spectralPhaseM2[i] += delta * delta2;
        }
    }

    void calculateGroupDelay() {
        double dwFactor = 0.5/(dF * twoPi<double>());
        groupDelay[0] = 2 * dwFactor * (spectralPhaseMean[1] - spectralPhaseMean[0]);
        for (int i = 1; i < Nfreq-1; i++) {
            groupDelay[i] = dwFactor * (spectralPhaseMean[i + 1] - spectralPhaseMean[i - 1]);
        }
        groupDelay[Nfreq - 1] = 2 * dwFactor * (spectralPhaseMean[Nfreq - 1] - spectralPhaseMean[Nfreq - 2]);
    }

public:
    spectralInterferometry() {
        setupFFT();
    }
    ~spectralInterferometry() {
        destroyFFT();
    }
    void resetFrequencies(size_t N, double fMin, double fMax) {
        if (N < 2) return;
        if (Nfreq == N && minF == fMin && maxF == fMax) return;
        frequencies = std::vector<double>(N);
        dF = (fMax - fMin) / (N - 1);
        Nfreq = N;
        minF = fMin;
        maxF = fMax;
        for (int i = 0; i < N; i++) {
            frequencies[i] = fMin + static_cast<double>(i) * dF;
        }
        if (FFTsize != Nfreq) {
            destroyFFT();
            setupFFT();
        }
        if (referenceDataA.size() > 0) {
            referenceDataAInterpolated = wavelengthToFrequency(frequencies, wavelengths, referenceDataA);
        }
        if (referenceDataB.size() > 0) {
            referenceDataAInterpolated = wavelengthToFrequency(frequencies, wavelengths, referenceDataB);
        }

        interferenceData = std::vector<double>(wavelengths);
        interferenceDataInterpolated = std::vector<double>(Nfreq, 0.0);
        spectralPhase = std::vector<double>(Nfreq, 0.0);
        spectralPhaseMean = spectralPhase;
        spectralPhaseM2 = spectralPhase;
        groupDelay = spectralPhase;
        phaseCount = 0;
        isConfigured = true;
    }

    void resetPhase() {
        phaseCount = 0;
        spectralPhaseMean = std::vector<double>(Nfreq, 0.0);
        spectralPhaseM2 = spectralPhaseMean;
    }

    void acquireNewPhase(OceanSpectrometer& s) {
        interferenceDataInterpolated = s.acquireSingleFrequency(frequencies);
        calculatePhase();
        updateWithNewPhase();
        calculateGroupDelay();
    }

    void acquireNewInterferogram(OceanSpectrometer& s) {
        interferenceDataInterpolated = s.acquireSingleFrequency(frequencies);
    }
    void calculatePhase() {
        for (int i = 0; i < Nfreq; i++) {
            hilbertRealBuffer[i] = interferenceDataInterpolated[i] - referenceDataBInterpolated[i] - referenceDataAInterpolated[i];
        }
        filteredHilbert(hilbertRealBuffer, hilbertRealBuffer, hilbertImagBuffer, filterT0, filterWidth, 8);
        spectralPhase = std::vector<double>(Nfreq);
        for (int i = 0; i < Nfreq; i++) {
            spectralPhase[i] = std::atan2(hilbertRealBuffer[i], hilbertImagBuffer[i]);
        }
        unwrap(spectralPhase);
    }
    double* getPhaseData() {
        return spectralPhase.data();
    }
    double* getMeanPhaseData() {
        return spectralPhaseMean.data();
    }
    double* getGroupDelayData() {
        return groupDelay.data();
    }
    double* getFrequencies() {
        return frequencies.data();
    }
    size_t getNfreq() {
        return Nfreq;
    }
    size_t getNtime() {
        return Nfreq / 2 + 1;
    }
    double* getInterferenceData() {
        return interferenceDataInterpolated.data();
    }
    double* getReferenceA() {
        return referenceDataAInterpolated.data();
    }
    double* getReferenceB() {
        return referenceDataBInterpolated.data();
    }
    std::vector<double>& getReferenceAVector() {
        return referenceDataAInterpolated;
    }
    std::vector<double>& getReferenceBVector() {
        return referenceDataBInterpolated;
    }

    double* getTimeData() {
        return fftDataReal.data();
    }

    double* getTimeReferenceA() {
        return fftReferenceAReal.data();
    }

    double* getTimeReferenceB() {
        return fftReferenceBReal.data();
    }
    double* getTimeScale() {
        return fftScale.data();
    }
    double* getTimeFilter() {
        return timeFilter.data();
    }
    void setTimeFilter(double t0, double sigma, double ord) {
        filterT0 = t0;
        filterWidth = sigma;
        filterOrder = ord;
    }
    void generateTimePlot() {
        fftw_execute_dft_r2c(fftPlanD2Z, interferenceDataInterpolated.data(), (fftw_complex*)hilbertTimeBuffer.data());
        fftw_execute_dft_r2c(fftPlanD2Z, referenceDataAInterpolated.data(), (fftw_complex*)fftReferenceA.data());
        fftw_execute_dft_r2c(fftPlanD2Z, referenceDataAInterpolated.data(), (fftw_complex*)fftReferenceB.data());
        double dt = 1.0e3 / (Nfreq * dF);
        double maxSignal = 0.0;
        for (int i = 0; i < (Nfreq / 2 + 1); i++) {

            fftReferenceAReal[i] = modSquared(fftReferenceA[i]);
            fftReferenceBReal[i] = modSquared(fftReferenceB[i]);
            fftDataReal[i] = modSquared(hilbertTimeBuffer[i] - fftReferenceA[i] - fftReferenceB[i]);
            maxSignal = maxN(fftDataReal[i], maxSignal);
            timeFilter[i] = std::exp(-std::pow(
                (static_cast<double>(i) * dt - filterT0) / filterWidth, filterOrder)
                / sqrt(2.0));
            fftScale[i] = static_cast<double>(i) * dt;
        }
        for (int i = 0; i < (Nfreq / 2 + 1); i++) {
            timeFilter[i] *= maxSignal;
        }
    }

    bool checkConfigurationStatus() {
        return isConfigured;
    }
    void acquireReferenceA(batchAcquisition& batchControl, const size_t N, const double integrationTime, const double secondsToWait, OceanSpectrometer& s) {
        if (N == 0) return;
        batchControl.acquireBatch(N, integrationTime, secondsToWait, s);
        wavelengths = s.wavelengthsCopy();
        referenceDataA = std::vector<double>(wavelengths.size(), 0.0);
        double normFactor = 1.0 / static_cast<double>(N);
        std::vector<double> fullSet = batchControl.getDataVector();
        for (int i = 0; i < wavelengths.size(); i++) {
            for (int j = 0; j < N; j++) {
                referenceDataA[i] += fullSet[i + wavelengths.size() * j];
            }
            referenceDataA[i] *= normFactor;
        }
        referenceDataAInterpolated = wavelengthToFrequency(frequencies, wavelengths, referenceDataA);
    }

    void acquireReferenceB(batchAcquisition& batchControl, const size_t N, const double integrationTime, const double secondsToWait, OceanSpectrometer& s) {
        if (N == 0) return;
        batchControl.acquireBatch(N, integrationTime, secondsToWait, s);
        wavelengths = s.wavelengthsCopy();
        referenceDataB = std::vector<double>(wavelengths.size(), 0.0);
        double normFactor = 1.0 / static_cast<double>(N);
        std::vector<double> fullSet = batchControl.getDataVector();
        for (int i = 0; i < wavelengths.size(); i++) {
            for (int j = 0; j < N; j++) {
                referenceDataB[i] +=fullSet[i + j *wavelengths.size()];
            }
            referenceDataB[i] *= normFactor;
        }
        referenceDataBInterpolated = wavelengthToFrequency(frequencies, wavelengths, referenceDataB);
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
        for (size_t j = 0; j < Nfreq; j++) {
            fs << frequencies[j];
            fs << " ";
            fs << spectralPhaseMean[j];
            fs << " ";
            fs << std::sqrt(spectralPhaseM2[j] / (phaseCount - 1));
            fs << '\x0A';
        }
    }
};
spectralInterferometry theInterferenceController;

//Main class for controlling the interface
class mainGui {
    bool queueUpdate = false;
    bool queueSliderUpdate = false;
    bool queueSliderMove = false;
    bool queueInterfaceValuesUpdate = false;
    bool isRunningLive = false;
    bool noSpectrometers = false;
    int sliderTarget = 0;
    
public:
    LweTextBox textBoxes[54];
    LweButton buttons[16];
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
    LweSpacer spacers[2];
    size_t pathTarget{};
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
        GtkCssProvider* textProvider = gtk_css_provider_new();
        gtk_css_provider_load_from_data(textProvider,
            "label, scale { font-family: Arial; font-weight: bold; }\n button, entry, textview { font-family: Arial; font-weight: bold; color: #EEEEEE; background-color: #191919; }", -1);
        gtk_style_context_add_provider_for_display(gdk_display_get_default(), GTK_STYLE_PROVIDER(textProvider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        GtkCssProvider* buttonShrinker = gtk_css_provider_new();
        gtk_css_provider_load_from_data(buttonShrinker,
            "label, scale, range, button, entry, textview { min-height: 10px; min-width: 10px; }", -1);
        gtk_style_context_add_provider_for_display(gdk_display_get_default(), GTK_STYLE_PROVIDER(buttonShrinker), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        buttons[0].init(("Run"), parentHandle, buttonCol1, 1, buttonWidth, 1, handleRunButton);
        buttons[1].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 3, buttonWidth/2, 1, handleGetOverlay0);
        buttons[2].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 4, buttonWidth / 2, 1, handleGetOverlay1);
        buttons[3].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 5, buttonWidth / 2, 1, handleGetOverlay2);
        buttons[3].init(("Acquire"), parentHandle, buttonCol1, 7, buttonWidth, 1, handleSave);

        buttons[4].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 3, buttonWidth / 2, 1, handleDeleteOverlay0);
        buttons[5].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 4, buttonWidth / 2, 1, handleDeleteOverlay1);
        buttons[6].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 5, buttonWidth / 2, 1, handleDeleteOverlay2);
        buttons[7].init(("\xf0\x9f\x95\xaf\xef\xb8\x8f"), parentHandle, buttonCol1, 2, buttonWidth / 2, 1, handleGetDarkSpectrum);
        buttons[8].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1 + buttonWidth / 2, 2, buttonWidth / 2, 1, handleDeleteDarkSpectrum);
        
        buttons[9].init(("Ref. A"), parentHandle, 0, 12, smallButton, 1, handleReferenceA);
        buttons[10].init(("Ref. B"), parentHandle, smallButton, 12, smallButton, 1, handleReferenceB);
        buttons[11].init(("Reset"), parentHandle, smallButton * 2, 12, smallButton, 1, handleResetPhase);
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
        buttons[7].init(("..."), parentHandle, buttonCol1 + buttonWidth / 2, 6, buttonWidth / 2, 1, saveFileDialogCallback, 0);
        textBoxes[48].init(window.parentHandle(4), 2, 0, 2, 1);
        textBoxes[49].init(window.parentHandle(4), 4, 0, 2, 1);
        textBoxes[50].init(window.parentHandle(4), 8, 0, 2, 1);
        textBoxes[51].init(window.parentHandle(4), 10, 0, 2, 1);
        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawSpectrum);
        checkBoxes[1].init(("Log"), window.parentHandle(4), 13, 0, 1, 1);
        buttons[10].init(("xlim"), window.parentHandle(4), 0, 0, 1, 1, handleRefreshRequest);
        buttons[10].setTooltip("Apply the entered x limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[10].squeeze();
        buttons[11].init(("ylim"), window.parentHandle(4), 6, 0, 1, 1, handleRefreshRequest);
        buttons[11].setTooltip("Apply the entered y limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[11].squeeze();
        buttons[12].init(("SVG"), window.parentHandle(3), 5, 0, 1, 1, svgCallback);
        buttons[12].setTooltip("Generate SVG files of the four line plots, with filenames based on the base path set above");
        buttons[12].squeeze();
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
        console.init(window.parentHandle(1), 0, 0, 1, 1);
        console.cPrint("Attached spectrometers:\n");
        
        g_signal_connect(window.window, "destroy", G_CALLBACK(destroyMainWindowCallback), NULL);
        window.present();
        initializeSpectrometers();
        pulldowns[0].init(parentHandle, 0, 0, 12, 1);
        drawBoxes[0].queueDraw();
        timeoutID = g_timeout_add(50, G_SOURCE_FUNC(updateDisplay), NULL);
    }

    void initializeSpectrometers() {
        int deviceCount = odapi_probe_devices();
        if (deviceCount == 0) {
            console.cPrint("I didn't find any spectrometers...\nMake sure you've installed the driver.\n");
            noSpectrometers = true;
            return;
        }
        int deviceIdCount = odapi_get_number_of_device_ids();

        std::vector<long> deviceIds(deviceIdCount);
        int retrievedIdCount = odapi_get_device_ids(deviceIds.data(), deviceIdCount);
        int error = 0;
        spectrometerSet = std::vector<OceanSpectrometer>(deviceIdCount);
        
		for (int i = 0; i < deviceIdCount; i++) {
			odapi_open_device(deviceIds[i], &error);
			// Get the device name
			const int nameLength = 32;
			char deviceName[nameLength] = { 0 };
			odapi_get_device_name(deviceIds[i], &error, deviceName, nameLength);
			if (error != 0) {
				console.cPrint("Failed to retrieve the spectrometer type. The error code is:  {}\n", error);
			}
			// and serial number
			int serialNumberLength = odapi_get_serial_number_maximum_length(deviceIds[i], &error);
			std::unique_ptr<char> serialNumber(new char[serialNumberLength]);
			odapi_get_serial_number(deviceIds[i], &error, serialNumber.get(), serialNumberLength);
			if (error != 0) {
				console.cPrint("Failed to retrieve the spectrometer serial number. The error code is:  {}\n", error);
			}
			else {
				console.cPrint("Device {}: {}\n    serial number: {}\n", i, deviceName, serialNumber.get());
                std::string newElement = Sformat("{}: {}", i, deviceName);
                pulldowns[0].addElement(newElement.c_str());
			}
            spectrometerSet[i].init(deviceIds[i]);
		}
	}

};
mainGui theGui;

void destroyMainWindowCallback() {
    for (int i = 0; i < spectrometerSet.size(); i++) {
        spectrometerSet[i].release();
    }
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

void handleGetOverlay0() {
    spectrometerSet[theGui.pulldowns[0].getValue()].acquireOverlay(0);
}

void handleGetOverlay1() {
    spectrometerSet[theGui.pulldowns[0].getValue()].acquireOverlay(1);
}

void handleGetOverlay2() {
    spectrometerSet[theGui.pulldowns[0].getValue()].acquireOverlay(2);
}

void handleDeleteOverlay0() {
    spectrometerSet[theGui.pulldowns[0].getValue()].deleteOverlay(0);
}

void handleDeleteOverlay1() {
    spectrometerSet[theGui.pulldowns[0].getValue()].deleteOverlay(1);
}

void handleDeleteOverlay2() {
    spectrometerSet[theGui.pulldowns[0].getValue()].deleteOverlay(2);
}

void handleGetDarkSpectrum() {
    spectrometerSet[theGui.pulldowns[0].getValue()].acquireDarkSpectrum();
}

void handleDeleteDarkSpectrum() {
    spectrometerSet[theGui.pulldowns[0].getValue()].disableDarkSpectrum();
}

void acquisitionThread(int activeSpectrometer, size_t N, double integrationTime, double secondsToWait, std::string path, bool timestamp) {
    theBatch.acquireBatch(N, integrationTime, secondsToWait, spectrometerSet[activeSpectrometer]);
    theBatch.save(path, timestamp);
    theGui.console.tPrint("Finished writing {}!\n", path);
}

void handleSave() {
    theGui.stopLive();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (spectrometerSet[activeSpectrometer].checkLock()) return;
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
    theInterferenceController.acquireReferenceA(theBatch, N, integrationTime, secondsToWait, spectrometerSet[activeSpectrometer]);
}

void referenceBAcquisitionThread(int activeSpectrometer, size_t N, double integrationTime, double secondsToWait) {
    theInterferenceController.acquireReferenceB(theBatch, N, integrationTime, secondsToWait, spectrometerSet[activeSpectrometer]);
}

void handleReferenceA() {
    handleResetController();
    if (!theInterferenceController.checkConfigurationStatus()) {
        theGui.console.cPrint("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stopLive();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (spectrometerSet[activeSpectrometer].checkLock()) return;
    size_t N = (size_t)theGui.textBoxes[14].valueDouble();
    double integrationTime = theGui.textBoxes[0].valueDouble();
    double waitTime = theGui.textBoxes[13].valueDouble();
    std::thread(referenceAAcquisitionThread, activeSpectrometer, N, integrationTime, waitTime).detach();
}

void handleReferenceB() {
    handleResetController();
    if (!theInterferenceController.checkConfigurationStatus()) {
        theGui.console.cPrint("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stopLive();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (spectrometerSet[activeSpectrometer].checkLock()) return;
    size_t N = (size_t)theGui.textBoxes[14].valueDouble();
    double integrationTime = theGui.textBoxes[0].valueDouble();
    double waitTime = theGui.textBoxes[13].valueDouble();
    std::thread(referenceBAcquisitionThread, activeSpectrometer, N, integrationTime, waitTime).detach();
}

void handleResetController() {
    size_t Nfreq = static_cast<size_t>(theGui.textBoxes[15].valueDouble());
    if (Nfreq < 2) return;
    double fMin = theGui.textBoxes[16].valueDouble();
    double fMax = theGui.textBoxes[17].valueDouble();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (fMax == 0.0) {
        fMax = spectrometerSet[activeSpectrometer].wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if (fMin == 0.0) {
        fMin = spectrometerSet[activeSpectrometer].wavelengths()[spectrometerSet[activeSpectrometer].size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    theInterferenceController.resetFrequencies(Nfreq, fMin, fMax);
}

void handleResetPhase() {
    theInterferenceController.resetPhase();
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
        drawInterferenceGroupDelay(area, cr, width, height, data);
        return;
    default:
        return;
    }
    if (theGui.noSpectrometersFound()) return;
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (!spectrometerSet[activeSpectrometer].initialized()) {
        theGui.console.cPrint("Not initialized - error {}\n",spectrometerSet[0].getErrorCode());
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
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_Spectrum.svg");
        sPlot.SVGPath = svgPath;
    }
    
    if (theGui.runningLive() && !spectrometerSet[activeSpectrometer].checkLock()) {
        spectrometerSet[activeSpectrometer].setIntegrationTime((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        spectrometerSet[activeSpectrometer].acquireSingle();
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
    sPlot.dataX = spectrometerSet[activeSpectrometer].wavelengths();
    sPlot.hasDataX = true;
    sPlot.data = spectrometerSet[activeSpectrometer].data();
    sPlot.Npts = spectrometerSet[activeSpectrometer].size();
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
    sPlot.xLabel = "Wavelength (nm)";
    sPlot.yLabel = "Spectrum (counts)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    sPlot.markers = false;
    if (spectrometerSet[activeSpectrometer].getOverlayCount() > 0) {
        int firstAdded = -1;
        int secondAdded = -1;
        sPlot.ExtraLines = spectrometerSet[activeSpectrometer].getOverlayCount();
        if (spectrometerSet[activeSpectrometer].hasOverlay0) {
            sPlot.data2 = spectrometerSet[activeSpectrometer].getOverlay(0);
            firstAdded = 0;
        }
        else if (spectrometerSet[activeSpectrometer].hasOverlay1) {
            sPlot.data2 = spectrometerSet[activeSpectrometer].getOverlay(1);
            firstAdded = 1;
        }
        else if (spectrometerSet[activeSpectrometer].hasOverlay2) {
            sPlot.data2 = spectrometerSet[activeSpectrometer].getOverlay(2);
            firstAdded = 2;
        }
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && spectrometerSet[activeSpectrometer].hasOverlay1) sPlot.data3 = spectrometerSet[activeSpectrometer].getOverlay(1);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && spectrometerSet[activeSpectrometer].hasOverlay2) sPlot.data3 = spectrometerSet[activeSpectrometer].getOverlay(2);
       
        if (sPlot.ExtraLines > 2) sPlot.data4 = spectrometerSet[activeSpectrometer].getOverlay(2);
    }
    sPlot.plot(cr);
}

void drawSpectrumFrequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {

    if (theGui.noSpectrometersFound()) return;
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (!spectrometerSet[activeSpectrometer].initialized()) {
        theGui.console.cPrint("Not initialized - error {}\n", spectrometerSet[0].getErrorCode());
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
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t Nfreq = static_cast<size_t>(theGui.textBoxes[15].valueDouble());
    if (Nfreq < 2) return;
    double fMin = theGui.textBoxes[16].valueDouble();
    double fMax = theGui.textBoxes[17].valueDouble();

    if (fMax == 0.0) {
        fMax = spectrometerSet[activeSpectrometer].wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if (fMin == 0.0) {
        fMin = spectrometerSet[activeSpectrometer].wavelengths()[spectrometerSet[activeSpectrometer].size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    double dF = (fMax - fMin) / static_cast<double>(Nfreq - 1);
    std::vector<double> frequencies(Nfreq);
    for (size_t i = 0; i < Nfreq; i++) {
        frequencies[i] = fMin + static_cast<double>(i) * dF;
    }
    std::vector<double> liveSpectrum;
    if (theGui.runningLive() && !spectrometerSet[activeSpectrometer].checkLock() && Nfreq > 0) {
        spectrometerSet[activeSpectrometer].setIntegrationTime((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        liveSpectrum = spectrometerSet[activeSpectrometer].acquireSingleFrequency(frequencies);
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
    sPlot.Npts = Nfreq;
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
    if (spectrometerSet[activeSpectrometer].getOverlayCount() > 0) {
        int firstAdded = -1;
        int secondAdded = -1;
        sPlot.ExtraLines = spectrometerSet[activeSpectrometer].getOverlayCount();
        if (spectrometerSet[activeSpectrometer].hasOverlay0) {
            sPlot.data2 = spectrometerSet[activeSpectrometer].getOverlayFrequency(0, frequencies);
            firstAdded = 0;
        }
        else if (spectrometerSet[activeSpectrometer].hasOverlay1) {
            sPlot.data2 = spectrometerSet[activeSpectrometer].getOverlayFrequency(1, frequencies);
            firstAdded = 1;
        }
        else if (spectrometerSet[activeSpectrometer].hasOverlay2) {
            sPlot.data2 = spectrometerSet[activeSpectrometer].getOverlayFrequency(2, frequencies);
            firstAdded = 2;
        }
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && spectrometerSet[activeSpectrometer].hasOverlay1) sPlot.data3 = spectrometerSet[activeSpectrometer].getOverlayFrequency(1, frequencies);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && spectrometerSet[activeSpectrometer].hasOverlay2) sPlot.data3 = spectrometerSet[activeSpectrometer].getOverlayFrequency(2, frequencies);

        if (sPlot.ExtraLines > 2) sPlot.data4 = spectrometerSet[activeSpectrometer].getOverlayFrequency(2, frequencies);
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
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[0].copyBuffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t Nfreq = static_cast<size_t>(theGui.textBoxes[15].valueDouble());
    if (Nfreq < 2) return;
    double fMin = theGui.textBoxes[16].valueDouble();
    double fMax = theGui.textBoxes[17].valueDouble();
    if (fMax == 0.0) {
        fMax = spectrometerSet[0].wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
        for (int i = 1; i < spectrometerSet.size(); i++) {
            fMax = maxN(fMax, 1e-3 * lightC<double>() / spectrometerSet[i].wavelengths()[0]);
        }
    }
    if (fMin == 0.0) {
        fMin = spectrometerSet[0].wavelengths()[spectrometerSet[0].size() - 1];
        fMin = 1e-3 * lightC<double>() / fMax;
        for (int i = 1; i < spectrometerSet.size(); i++) {
            fMin = minN(fMin, 1e-3 * lightC<double>() / spectrometerSet[i].wavelengths()[spectrometerSet[0].size() - 1]);
        }
    }

    double dF = (fMax - fMin) / static_cast<double>(Nfreq - 1);
    std::vector<double> frequencies(Nfreq);
    for (size_t i = 0; i < Nfreq; i++) {
        frequencies[i] = fMin + static_cast<double>(i) * dF;
    }
    std::vector<std::vector<double>> liveSpectra(spectrometerSet.size());
    if (theGui.runningLive() && Nfreq > 0) {
        for (int i = 0; i < spectrometerSet.size(); i++) {
            liveSpectra[i] = spectrometerSet[i].acquireSingleFrequency(frequencies);
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
    sPlot.Npts = Nfreq;
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
    sPlot.ExtraLines = liveSpectra.size()-1;
    if (liveSpectra.size() > 0) sPlot.data2 = liveSpectra[1].data();
    if (liveSpectra.size() > 1) sPlot.data3 = liveSpectra[2].data();
    if (liveSpectra.size() > 2) sPlot.data4 = liveSpectra[3].data();
    sPlot.plot(cr);
}

void drawInterferenceSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    handleResetController();
    if (!theInterferenceController.checkConfigurationStatus()) handleResetController();
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

    if (theGui.runningLive() && !spectrometerSet[theGui.pulldowns[0].getValue()].checkLock()) {
        spectrometerSet[theGui.pulldowns[0].getValue()].setIntegrationTime((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.acquireNewInterferogram(spectrometerSet[theGui.pulldowns[0].getValue()]);
    }
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.getFrequencies();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.getInterferenceData();
    sPlot.Npts = theInterferenceController.getNfreq();
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
    if (theInterferenceController.getReferenceAVector().size() == theInterferenceController.getNfreq()) {
        sPlot.data2 = theInterferenceController.getReferenceA();
        sPlot.ExtraLines = 1;
        if (theInterferenceController.getReferenceBVector().size() == theInterferenceController.getNfreq()) {
            sPlot.data3 = theInterferenceController.getReferenceB();
            sPlot.ExtraLines = 2;
        }
    }
    sPlot.plot(cr);
}

void drawInterferenceSpectrumTime(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
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
    if (theGui.runningLive() && !spectrometerSet[theGui.pulldowns[0].getValue()].checkLock()) {
        spectrometerSet[theGui.pulldowns[0].getValue()].setIntegrationTime((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.acquireNewInterferogram(spectrometerSet[theGui.pulldowns[0].getValue()]);
    }
    theInterferenceController.setTimeFilter(t0, sigma, ord);
    theInterferenceController.generateTimePlot();

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.getTimeScale();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.getTimeData();
    sPlot.Npts = theInterferenceController.getNtime();
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
    sPlot.data2 = theInterferenceController.getTimeReferenceA();
    sPlot.data3 = theInterferenceController.getTimeReferenceB();
    sPlot.data4 = theInterferenceController.getTimeFilter();
    sPlot.ExtraLines = 3;
    sPlot.plot(cr);
}


void drawInterferencePhase(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
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

    if (theGui.runningLive() && !spectrometerSet[theGui.pulldowns[0].getValue()].checkLock()) {
        spectrometerSet[theGui.pulldowns[0].getValue()].setIntegrationTime((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.acquireNewPhase(spectrometerSet[theGui.pulldowns[0].getValue()]);
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.getFrequencies();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.getMeanPhaseData();
    sPlot.Npts = theInterferenceController.getNfreq();
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

void drawInterferenceGroupDelay(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
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

    if (theGui.runningLive() && !spectrometerSet[theGui.pulldowns[0].getValue()].checkLock()) {
        spectrometerSet[theGui.pulldowns[0].getValue()].setIntegrationTime((unsigned long)round(1000 * theGui.textBoxes[0].valueDouble()));
        theInterferenceController.acquireNewPhase(spectrometerSet[theGui.pulldowns[0].getValue()]);
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = theInterferenceController.getFrequencies();
    sPlot.hasDataX = true;
    sPlot.data = theInterferenceController.getGroupDelayData();
    sPlot.Npts = theInterferenceController.getNfreq();
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
void pathFromDialogBox(GtkDialog* dialog, int response) {
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        std::string s(g_file_get_path(file));
        theGui.filePaths[theGui.pathTarget].overwritePrint("{}", s);
    }
    g_object_unref(dialog);
    }
void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget) {
    theGui.pathTarget = (size_t)pathTarget;
    //get around bug in GTK4 by opening dialog box directly in cocoa on mac
#ifdef __APPLE__
    NSString* filePath;
    NSSavePanel* savePanel = [NSSavePanel savePanel];
    if ([savePanel runModal] == NSModalResponseOK) {
        filePath = [savePanel URL].path;
        theGui.filePaths[theGui.pathTarget].overwritePrint("{}", [filePath UTF8String]);
    }
    return;
#else
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new("Save File", theGui.window.windowHandle(), GTK_FILE_CHOOSER_ACTION_SAVE, "Ok", "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
#endif
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
    GtkApplication* app = gtk_application_new("nickkarpowicz.lightwave", (GApplicationFlags)0);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    return g_application_run(G_APPLICATION(app), argc, argv);
}
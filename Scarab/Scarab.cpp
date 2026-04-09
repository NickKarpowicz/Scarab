#include <thread>

#include "api/OceanDirectAPI.h"
#include "ExternalLibraries/LightwaveExplorerGraphicalClasses.h"
#include "spectrometer.hpp"
#include "ocean_spectrometer.hpp"
#include "batch_acquisition.hpp"
#include "interferometry.hpp"

bool update_display();
void handle_run_button();
void draw_spectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void draw_spectrum_frequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void draw_spectra_frequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void draw_interference_spectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void draw_interference_spectrum_time(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void draw_interference_phase(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void draw_interference_group_delay(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drop_down_change_callback();
void handle_get_overlay0();
void handle_get_overlay1();
void handle_get_overlay2();
void handle_collapse_panel();
void handle_delete_overlay0();
void handle_delete_overlay1();
void handle_delete_overlay2();
void handle_get_dark_spectrum();
void handle_delete_dark_spectrum();
void handle_refresh_request();
void save_path_callback();
void handle_save();
void handle_reference_A();
void handle_reference_B();
void handle_reset_controller();
void handle_reset_phase();
void handle_save_phase();
void svg_callback();

std::vector<std::unique_ptr<Spectrometer>> spectrometerSet;
BatchAcquisition theBatch;
SpectralInterferometry the_interference_controller;

//Main class for controlling the interface
class MainGui {
    bool queue_update = false;
    bool queue_interface_values_update = false;
    bool queue_save_path_update = false;
    bool is_running_live = false;
public:
    LweTextBox text_boxes[54];
    LweButton buttons[18];
    LweConsole console;
    LweConsole sequence;
    LweTextBox file_paths[4];
    LwePulldown pulldowns[10];
    LweDrawBox draw_boxes[8];
    LweCheckBox checkboxes[4];
    LweWindow window;
    std::string path_buffer;
    int saveSVG = 0;
    bool loadedDefaults = false;
    unsigned int timeoutID = 0;
    void requestPlotUpdate() {
        queue_update = true;
    }

    void request_live() {
        is_running_live = true;
    }

    void stop_live() {
        is_running_live = false;
    }

    [[nodiscard]] constexpr bool running_live() { return is_running_live; }
    void apply_update() {
        if (queue_update) {
            queue_update = false;
            draw_boxes[0].queueDraw();
        }
        if (queue_save_path_update && path_buffer != "?LWE_LOADING??") {
            queue_save_path_update = false;
            if (path_buffer == "?LWE_NOPATH??") return;
            file_paths[0].overwritePrint(path_buffer);
        }
    }


    void requestInterfaceValuesUpdate() {
        queue_interface_values_update = true;
    }
    void requestSavePathUpdate() {
#ifndef __APPLE__
        path_buffer = std::string("?LWE_LOADING??");
        pathFromSaveDialog(path_buffer, "txt", "Text (.txt)");
#endif
        queue_save_path_update = true;
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

        buttons[0].init(("Run"), parentHandle, buttonCol1, 1, buttonWidth, 1, handle_run_button);
        buttons[1].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 3, buttonWidth/2, 1, handle_get_overlay0);
        buttons[2].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 4, buttonWidth / 2, 1, handle_get_overlay1);
        buttons[3].init(("\xf0\x9f\x93\x88"), parentHandle, buttonCol1, 5, buttonWidth / 2, 1, handle_get_overlay2);
        buttons[3].init(("Acquire"), parentHandle, buttonCol1, 7, buttonWidth, 1, handle_save);

        buttons[4].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 3, buttonWidth / 2, 1, handle_delete_overlay0);
        buttons[5].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 4, buttonWidth / 2, 1, handle_delete_overlay1);
        buttons[6].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 5, buttonWidth / 2, 1, handle_delete_overlay2);
        buttons[7].init(("\xf0\x9f\x95\xaf\xef\xb8\x8f"), parentHandle, buttonCol1, 2, buttonWidth / 2, 1, handle_get_dark_spectrum);
        buttons[8].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1 + buttonWidth / 2, 2, buttonWidth / 2, 1, handle_delete_dark_spectrum);

        buttons[9].init(("Ref. A"), parentHandle, 0, 12, smallButton, 1, handle_reference_A);
        buttons[10].init(("Ref. B"), parentHandle, smallButton, 12, smallButton, 1, handle_reference_B);
        buttons[11].init(("Reset"), parentHandle, smallButton * 2, 12, smallButton, 1, handle_reset_phase);
        buttons[12].init(("Save"), parentHandle, smallButton * 3, 12, smallButton, 1, handle_save_phase);

        //RGB active
        text_boxes[1].init(parentHandle, textCol1, 2, textWidth, 1);
        text_boxes[2].init(parentHandle, textCol2, 2, textWidth, 1);
        text_boxes[3].init(parentHandle, textCol3, 2, textWidth, 1);

        //RGB overlay0
        text_boxes[4].init(parentHandle, textCol1, 3, textWidth, 1);
        text_boxes[5].init(parentHandle, textCol2, 3, textWidth, 1);
        text_boxes[6].init(parentHandle, textCol3, 3, textWidth, 1);

        //RGB overlay1
        text_boxes[7].init(parentHandle, textCol1, 4, textWidth, 1);
        text_boxes[8].init(parentHandle, textCol2, 4, textWidth, 1);
        text_boxes[9].init(parentHandle, textCol3, 4, textWidth, 1);

        //RGB overlay2
        text_boxes[10].init(parentHandle, textCol1, 5, textWidth, 1);
        text_boxes[11].init(parentHandle, textCol2, 5, textWidth, 1);
        text_boxes[12].init(parentHandle, textCol3, 5, textWidth, 1);

        text_boxes[0].init(parentHandle, textCol3, 7, textWidth, 1);
        text_boxes[0].setLabel(-3 * textWidth, 0, "Exposure (ms)");
        text_boxes[0].overwritePrint(std::string("40"));

        text_boxes[13].init(parentHandle, textCol3, 8, textWidth, 1);
        text_boxes[13].setLabel(-3 * textWidth, 0, "Pause, count (s, #)");
        text_boxes[14].init(parentHandle, textCol4, 8, textWidth, 1);
        text_boxes[13].overwritePrint(std::string("0"));
        text_boxes[14].overwritePrint(std::string("10"));


        text_boxes[15].init(parentHandle, textCol3, 9, 2, 1);
        text_boxes[15].setMaxCharacters(6);
        text_boxes[16].init(parentHandle, textCol4, 9, 2, 1);
        text_boxes[16].setMaxCharacters(6);
        text_boxes[17].init(parentHandle, textCol5, 9, 2, 1);
        text_boxes[17].setMaxCharacters(6);
        text_boxes[15].setLabel(-3 * textWidth, 0, "Freqs. (#, min,max)");
        text_boxes[15].overwritePrint(std::string("2048"));

        text_boxes[18].init(parentHandle, textCol3, 11, 2, 1);
        text_boxes[18].setMaxCharacters(6);
        text_boxes[19].init(parentHandle, textCol4, 11, 2, 1);
        text_boxes[19].setMaxCharacters(6);
        text_boxes[20].init(parentHandle, textCol5, 11, 2, 1);
        text_boxes[20].setMaxCharacters(6);
        text_boxes[18].setLabel(-3 * textWidth, 0, "Filter (t0, sig., ord.)");
        text_boxes[18].overwritePrint(std::string("300"));
        text_boxes[19].overwritePrint(std::string("80"));
        text_boxes[20].overwritePrint(std::string("12"));

        file_paths[0].init(parentHandle, 0, 6, 10, 1);
        file_paths[0].setMaxCharacters(pathChars);
        file_paths[0].overwritePrint(Sformat("DefaultOutput.txt"));
        checkboxes[2].init("\xe2\x8c\x9a", parentHandle, buttonCol1 + buttonWidth/2, 8, 2, 1);
        buttons[7].init(("..."), parentHandle, buttonCol1 + buttonWidth / 2, 6, buttonWidth / 2, 1, save_path_callback, 0);
        text_boxes[48].init(window.parentHandle(4), 3, 0, 2, 1);
        text_boxes[49].init(window.parentHandle(4), 5, 0, 2, 1);
        text_boxes[50].init(window.parentHandle(4), 9, 0, 2, 1);
        text_boxes[51].init(window.parentHandle(4), 11, 0, 2, 1);
        draw_boxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        draw_boxes[0].setDrawingFunction(draw_spectrum);
        checkboxes[1].init(("Log"), window.parentHandle(4), 13, 0, 1, 1);
        checkboxes[3].init(("Average phase"), window.parentHandle(4), 15, 0, 1, 1);
        buttons[13].init(("xlim"), window.parentHandle(4), 2, 0, 1, 1, handle_refresh_request);
        buttons[13].setTooltip("Apply the entered x limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[13].squeeze();
        buttons[14].init(("ylim"), window.parentHandle(4), 7, 0, 1, 1, handle_refresh_request);
        buttons[14].setTooltip("Apply the entered y limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[14].squeeze();
        buttons[15].init(("SVG"), window.parentHandle(4), 1, 0, 1, 1, svg_callback);
        buttons[15].setTooltip("Generate SVG files of the four line plots, with filenames based on the base path set above");
        buttons[15].squeeze();
        buttons[16].init("\xe2\x86\x94\xef\xb8\x8f", window.parentHandle(4), 0, 0, 1, 1, handle_collapse_panel);
        buttons[16].setTooltip("Collapse/expand the data entry panel");

        for (int i = 0; i < 18; i++) {
            text_boxes[i].setMaxCharacters(6);
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
        gtk_widget_set_visible(text_boxes[18].label, false);
        gtk_widget_set_visible(text_boxes[18].elementHandle, false);
        gtk_widget_set_visible(text_boxes[19].elementHandle, false);
        gtk_widget_set_visible(text_boxes[20].elementHandle, false);
        console.init(window.parentHandle(1), 0, 0, 1, 1);
        console.cPrint("Attached spectrometers:\n");

        window.present();
        initializeSpectrometers();
        pulldowns[0].init(parentHandle, 0, 0, 12, 1);
        draw_boxes[0].queueDraw();
        timeoutID = g_timeout_add(20, G_SOURCE_FUNC(update_display), NULL);
    }

    void initializeSpectrometers() {
        int deviceCount = odapi_probe_devices();
        if (deviceCount == 0) {
            console.cPrint("I didn't find any spectrometers...\nMake sure you've installed the driver.\n");
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
                return;
            }
			// Get the device name
			const int nameLength = 128;
			char deviceName[nameLength] = { 0 };
			odapi_get_device_name(device_ids[i], &error, deviceName, nameLength);
			if (error != 0) {
				console.cPrint("Failed to retrieve the spectrometer type. The error code is:  {}\n", error);
                return;
			}
			// and serial number
			int serialNumberLength = odapi_get_serial_number_maximum_length(device_ids[i], &error);
            if (error != 0) {
                console.cPrint("Failed to retrieve the serial number length. The error code is:  {}\n", error);
                return;
            }
            std::unique_ptr<char> serialNumber(new char[serialNumberLength + 1] {});
			odapi_get_serial_number(device_ids[i], &error, serialNumber.get(), serialNumberLength);
			if (error != 0) {
				console.cPrint("Failed to retrieve the spectrometer serial number. The error code is:  {}\n", error);
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



void handle_run_button() {
    if (theGui.running_live()) {
        theGui.stop_live();
    }
    else {
        theGui.request_live();
    }
}

void svg_callback() {
    theGui.saveSVG = 1;
    theGui.requestPlotUpdate();
}

void handle_refresh_request() {
    theGui.requestPlotUpdate();
}

void handle_get_overlay0() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_overlay(0);
}

void drop_down_change_callback(){
    bool visibility = theGui.pulldowns[1].getValue()>2;
    gtk_widget_set_visible(theGui.buttons[9].elementHandle, visibility);
    gtk_widget_set_visible(theGui.buttons[10].elementHandle, visibility);
    gtk_widget_set_visible(theGui.buttons[11].elementHandle, visibility);
    gtk_widget_set_visible(theGui.buttons[12].elementHandle, visibility);
    gtk_widget_set_visible(theGui.text_boxes[18].label, visibility);
    gtk_widget_set_visible(theGui.text_boxes[18].elementHandle, visibility);
    gtk_widget_set_visible(theGui.text_boxes[19].elementHandle, visibility);
    gtk_widget_set_visible(theGui.text_boxes[20].elementHandle, visibility);
}

void handle_get_overlay1() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_overlay(1);
}

void handle_get_overlay2() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_overlay(2);
}

void handle_collapse_panel() {
    theGui.window.toggleSettingsPanel();
}

void handle_delete_overlay0() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).delete_overlay(0);
}

void handle_delete_overlay1() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).delete_overlay(1);
}

void handle_delete_overlay2() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).delete_overlay(2);
}

void handle_get_dark_spectrum() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).acquire_dark_spectrum();
}

void handle_delete_dark_spectrum() {
    (*spectrometerSet[theGui.pulldowns[0].getValue()]).disabledark_spectrum();
}

void acquisitionThread(int activeSpectrometer, size_t N, double integration_time, double seconds_to_wait, std::string path, bool timestamp) {
    theBatch.acquire_batch(N, integration_time, seconds_to_wait, (*spectrometerSet[activeSpectrometer]));
    theBatch.save(path, timestamp);
    theGui.console.tPrint("Finished writing {}!\n", path);
}

void handle_save() {
    theGui.stop_live();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((*spectrometerSet[activeSpectrometer]).checkLock()) return;
    std::string path;
    theGui.file_paths[0].copyBuffer(path);
    bool timestamp = theGui.checkboxes[2].isChecked();
    size_t N = (size_t)theGui.text_boxes[14].valueDouble();
    double integration_time = theGui.text_boxes[0].valueDouble();
    double waitTime = theGui.text_boxes[13].valueDouble();
    std::thread(acquisitionThread, activeSpectrometer, N, integration_time, waitTime, path, timestamp).detach();
}

void handle_save_phase() {
    std::string path;
    theGui.file_paths[0].copyBuffer(path);
    bool timestamp = theGui.checkboxes[2].isChecked();
    the_interference_controller.save(path, timestamp);
}

void referenceAAcquisitionThread(int activeSpectrometer, size_t N, double integration_time, double seconds_to_wait) {
    the_interference_controller.acquire_reference_A(theBatch, N, integration_time, seconds_to_wait, (*spectrometerSet[activeSpectrometer]));
}

void referenceBAcquisitionThread(int activeSpectrometer, size_t N, double integration_time, double seconds_to_wait) {
    the_interference_controller.acquireReferenceB(theBatch, N, integration_time, seconds_to_wait, (*spectrometerSet[activeSpectrometer]));
}

void handle_reference_A() {
    handle_reset_controller();
    if (!the_interference_controller.check_configuration_status()) {
        theGui.console.cPrint("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stop_live();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((*spectrometerSet[activeSpectrometer]).checkLock()) return;
    size_t N = (size_t)theGui.text_boxes[14].valueDouble();
    double integration_time = theGui.text_boxes[0].valueDouble();
    double waitTime = theGui.text_boxes[13].valueDouble();
    std::thread(referenceAAcquisitionThread, activeSpectrometer, N, integration_time, waitTime).detach();
}

void handle_reference_B() {
    handle_reset_controller();
    if (!the_interference_controller.check_configuration_status()) {
        theGui.console.cPrint("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stop_live();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((*spectrometerSet[activeSpectrometer]).checkLock()) return;
    size_t N = (size_t)theGui.text_boxes[14].valueDouble();
    double integration_time = theGui.text_boxes[0].valueDouble();
    double waitTime = theGui.text_boxes[13].valueDouble();
    std::thread(referenceBAcquisitionThread, activeSpectrometer, N, integration_time, waitTime).detach();
}

void handle_reset_controller() {
    size_t num_freqs = static_cast<size_t>(theGui.text_boxes[15].valueDouble());
    if (num_freqs < 2) return;
    double fMin = theGui.text_boxes[16].valueDouble();
    double fMax = theGui.text_boxes[17].valueDouble();
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if (fMax == 0.0) {
        fMax = (*spectrometerSet[activeSpectrometer]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if (fMin == 0.0) {
        fMin = (*spectrometerSet[activeSpectrometer]).wavelengths()[(*spectrometerSet[activeSpectrometer]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    the_interference_controller.reset_frequencies(num_freqs, fMin, fMax);
}

void handle_reset_phase() {
    the_interference_controller.reset_phase();
}

void draw_spectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    switch (theGui.pulldowns[1].getValue()) {
    case 0:
        break;
    case 1:
        draw_spectrum_frequency(area, cr, width, height, data);
        return;
    case 2:
        draw_spectra_frequency(area, cr, width, height, data);
        return;
    case 3:
        draw_interference_spectrum(area, cr, width, height, data);
        return;
    case 4:
        draw_interference_spectrum_time(area, cr, width, height, data);
        return;
    case 5:
        draw_interference_phase(area, cr, width, height, data);
        return;
    case 6:
        draw_interference_group_delay(area, cr, width, height, data);
        return;
    default:
        return;
    }

    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).initialized()) {
        theGui.console.cPrint("Not initialized - error {}\n",(*spectrometerSet[0]).get_error_code());
        return;
    }

    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;

    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
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
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_Spectrum.svg");
        sPlot.SVGPath = svgPath;
    }

    if (theGui.running_live() && (activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).checkLock()) {
        (*spectrometerSet[activeSpectrometer]).set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].valueDouble()));
        (*spectrometerSet[activeSpectrometer]).acquire_single();
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
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
        if ((*spectrometerSet[activeSpectrometer]).has_overlay_0) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay(0);
            firstAdded = 0;
        }
        else if ((*spectrometerSet[activeSpectrometer]).has_overlay_1) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay(1);
            firstAdded = 1;
        }
        else if ((*spectrometerSet[activeSpectrometer]).has_overlay_2) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay(2);
            firstAdded = 2;
        }
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && (*spectrometerSet[activeSpectrometer]).has_overlay_1) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlay(1);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && (*spectrometerSet[activeSpectrometer]).has_overlay_2) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlay(2);

        if (sPlot.ExtraLines > 2) sPlot.data4 = (*spectrometerSet[activeSpectrometer]).get_overlay(2);
    }
    sPlot.plot(cr);
}

void draw_spectrum_frequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    int activeSpectrometer = theGui.pulldowns[0].getValue();
    if ((activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).initialized()) {
        theGui.console.cPrint("Not initialized - error {}\n", (*spectrometerSet[0]).get_error_code());
        return;
    }

    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;

    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
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
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t num_freqs = static_cast<size_t>(theGui.text_boxes[15].valueDouble());
    if (num_freqs < 2) return;
    double fMin = theGui.text_boxes[16].valueDouble();
    double fMax = theGui.text_boxes[17].valueDouble();

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
    if (theGui.running_live() && (activeSpectrometer < spectrometerSet.size()) && !(*spectrometerSet[activeSpectrometer]).checkLock() && num_freqs > 0) {
        (*spectrometerSet[activeSpectrometer]).set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].valueDouble()));
        liveSpectrum = (*spectrometerSet[activeSpectrometer]).acquire_single_frequency(frequencies);
    }
    else return;

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
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
        if ((*spectrometerSet[activeSpectrometer]).has_overlay_0) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay_frequency(0, frequencies);
            firstAdded = 0;
        }
        else if ((*spectrometerSet[activeSpectrometer]).has_overlay_1) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay_frequency(1, frequencies);
            firstAdded = 1;
        }
        else if ((*spectrometerSet[activeSpectrometer]).has_overlay_2) {
            sPlot.data2 = (*spectrometerSet[activeSpectrometer]).get_overlay_frequency(2, frequencies);
            firstAdded = 2;
        }
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && (*spectrometerSet[activeSpectrometer]).has_overlay_1) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlay_frequency(1, frequencies);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && (*spectrometerSet[activeSpectrometer]).has_overlay_2) sPlot.data3 = (*spectrometerSet[activeSpectrometer]).get_overlay_frequency(2, frequencies);

        if (sPlot.ExtraLines > 2) sPlot.data4 = (*spectrometerSet[activeSpectrometer]).get_overlay_frequency(2, frequencies);
    }
    sPlot.plot(cr);
}

void draw_spectra_frequency(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;

    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
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
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t num_freqs = static_cast<size_t>(theGui.text_boxes[15].valueDouble());
    if (num_freqs < 2) return;
    double fMin = theGui.text_boxes[16].valueDouble();
    double fMax = theGui.text_boxes[17].valueDouble();
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
    if (theGui.running_live() && num_freqs > 0) {
        for (int i = 0; i < spectrometerSet.size(); i++) {
            liveSpectra[i] = (*spectrometerSet[i]).acquire_single_frequency(frequencies);
        }
    }
    else return;

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
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

void draw_interference_spectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    handle_reset_controller();
    if (!the_interference_controller.check_configuration_status()) handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_InterferenceF.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
    }

    if (theGui.running_live() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].valueDouble()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].isChecked());
        the_interference_controller.acquire_new_interferogram((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = the_interference_controller.get_frequencies();
    sPlot.hasDataX = true;
    sPlot.data = the_interference_controller.get_interference_data();
    sPlot.Npts = the_interference_controller.get_num_freqs();
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
    if (the_interference_controller.get_reference_A_vector().size() == the_interference_controller.get_num_freqs()) {
        sPlot.data2 = the_interference_controller.get_reference_A();
        sPlot.ExtraLines = 1;
        if (the_interference_controller.get_reference_B_vector().size() == the_interference_controller.get_num_freqs()) {
            sPlot.data3 = the_interference_controller.get_reference_B();
            sPlot.ExtraLines = 2;
        }
    }
    sPlot.plot(cr);
}

void draw_interference_spectrum_time(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_InterferenceF.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
    }
    double t0 = theGui.text_boxes[18].valueDouble();
    double sigma = theGui.text_boxes[19].valueDouble();
    double ord = theGui.text_boxes[20].valueDouble();
    if (theGui.running_live() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].valueDouble()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].isChecked());
        the_interference_controller.acquire_new_interferogram((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }
    the_interference_controller.set_time_filter(t0, sigma, ord);
    the_interference_controller.generate_time_plot();

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = the_interference_controller.get_time_scale();
    sPlot.hasDataX = true;
    sPlot.data = the_interference_controller.get_time_data();
    sPlot.Npts = the_interference_controller.get_num_times();
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
    sPlot.data2 = the_interference_controller.get_time_reference_A();
    sPlot.data3 = the_interference_controller.get_time_reference_B();
    sPlot.data4 = the_interference_controller.get_time_filter();
    sPlot.ExtraLines = 3;
    sPlot.plot(cr);
}


void draw_interference_phase(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_InterferencePhase.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
    }

    if (theGui.running_live() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].valueDouble()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].isChecked());
        the_interference_controller.acquire_new_phase((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = the_interference_controller.get_frequencies();
    sPlot.hasDataX = true;
    sPlot.data = the_interference_controller.get_mean_phase_data();
    sPlot.Npts = the_interference_controller.get_num_freqs();
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

void draw_interference_group_delay(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkboxes[1].isChecked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[48].valueDouble();
    double xMax = theGui.text_boxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[50].valueDouble();
    double yMax = theGui.text_boxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copyBuffer(svgPath);
        svgPath.append("_InterferencePhase.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[1].valueDouble() + theGui.text_boxes[2].valueDouble() + theGui.text_boxes[3].valueDouble())) {
        mainColor = LweColor(theGui.text_boxes[1].valueDouble(), theGui.text_boxes[2].valueDouble(), theGui.text_boxes[3].valueDouble(), 1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if (0.0 != (theGui.text_boxes[4].valueDouble() + theGui.text_boxes[5].valueDouble() + theGui.text_boxes[6].valueDouble())) {
        overLay0Color = LweColor(theGui.text_boxes[4].valueDouble(), theGui.text_boxes[5].valueDouble(), theGui.text_boxes[6].valueDouble(), 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if (0.0 != (theGui.text_boxes[7].valueDouble() + theGui.text_boxes[8].valueDouble() + theGui.text_boxes[9].valueDouble())) {
        overLay1Color = LweColor(theGui.text_boxes[7].valueDouble(), theGui.text_boxes[8].valueDouble(), theGui.text_boxes[9].valueDouble(), 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if (0.0 != (theGui.text_boxes[10].valueDouble() + theGui.text_boxes[11].valueDouble() + theGui.text_boxes[12].valueDouble())) {
        overLay2Color = LweColor(theGui.text_boxes[10].valueDouble(), theGui.text_boxes[11].valueDouble(), theGui.text_boxes[12].valueDouble(), 1);
    }

    if (theGui.running_live() && ! (*spectrometerSet[theGui.pulldowns[0].getValue()]).checkLock()) {
        (*spectrometerSet[theGui.pulldowns[0].getValue()]).set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].valueDouble()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].isChecked());
        the_interference_controller.acquire_new_phase((*spectrometerSet[theGui.pulldowns[0].getValue()]));
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = the_interference_controller.get_frequencies();
    sPlot.hasDataX = true;
    sPlot.data = the_interference_controller.get_group_delay_data();
    sPlot.Npts = the_interference_controller.get_num_freqs();
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

bool update_display() {
    theGui.requestPlotUpdate();
    theGui.console.updateFromBuffer();
    theGui.apply_update();
    return true;
}

void save_path_callback() {
#ifdef __APPLE__
    theGui.path_buffer = pathFromAppleSaveDialog();
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

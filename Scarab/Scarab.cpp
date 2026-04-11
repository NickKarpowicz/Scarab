#include "Scarab.hpp"
#include "batch_acquisition.hpp"
#include "gtk_wrappers.hpp"
#include "interferometry.hpp"
#include "ocean_spectrometer.hpp"
#include "spectrometer.hpp"
#include <thread>

std::vector<std::unique_ptr<Spectrometer>> spectrometer_set;
BatchAcquisition theBatch;
SpectralInterferometry the_interference_controller;

// Main class for controlling the interface
class MainGui {
    bool queue_update = false;
    bool queue_interface_values_update = false;
    bool queue_save_path_update = false;
    bool is_running_live = false;

  public:
    GtkTextBox text_boxes[26];
    GtkPushbutton buttons[19];
    GtkConsole console;
    GtkTextBox file_paths[2];
    GtkPulldown pulldowns[4];
    GtkDrawBox draw_boxes[4];
    GtkCheck checkboxes[4];
    GtkMainWindow window;
    std::string path_buffer;
    int saveSVG = 0;
    bool loadedDefaults = false;
    unsigned int timeoutID = 0;
    void requestPlotUpdate() { queue_update = true; }

    void request_live() { is_running_live = true; }

    void stop_live() { is_running_live = false; }

    [[nodiscard]] constexpr bool running_live() { return is_running_live; }
    void apply_update() {
        if(queue_update) {
            queue_update = false;
            draw_boxes[0].queue_draw();
        }
        if(queue_save_path_update && path_buffer != "?LWE_LOADING??") {
            queue_save_path_update = false;
            if(path_buffer == "?LWE_NOPATH??")
                return;
            file_paths[0].overwrite_print(path_buffer);
        }
    }

    void requestInterfaceValuesUpdate() { queue_interface_values_update = true; }
    void requestSavePathUpdate() {
#ifndef __APPLE__
        path_buffer = std::string("?LWE_LOADING??");
        path_from_save_dialog(path_buffer, "txt", "Text (.txt)");
#endif
        queue_save_path_update = true;
    }
    void activate(GtkApplication *app) {
        int buttonWidth = 4;
        int smallButton = 3;
        int textWidth = 2;
        int plotWidth = 12;
        int plotHeight = 6;
        int pathChars = 42;
        int textCol1 = textWidth;
        int textCol2 = 2 * textWidth;
        int textCol3 = 3 * textWidth;
        int textCol4 = 4 * textWidth;
        int textCol5 = 5 * textWidth;
        int buttonCol1 = 8;

        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", true, NULL);
        window.init(app, "Spectrometer Control and Recording for Attosecond Beamlines", 1400, 800);
        GtkWidget *parent_handle = window.parent_handle();

        buttons[0].init(("Run"), parent_handle, buttonCol1, 1, buttonWidth, 1, handle_run_button);
        buttons[1].init(("\xf0\x9f\x93\x88"),
                        parent_handle,
                        buttonCol1,
                        3,
                        buttonWidth / 2,
                        1,
                        handle_get_overlay0);
        buttons[2].init(("\xf0\x9f\x93\x88"),
                        parent_handle,
                        buttonCol1,
                        4,
                        buttonWidth / 2,
                        1,
                        handle_get_overlay1);
        buttons[3].init(("\xf0\x9f\x93\x88"),
                        parent_handle,
                        buttonCol1,
                        5,
                        buttonWidth / 2,
                        1,
                        handle_get_overlay2);
        buttons[18].init(("Acquire"), parent_handle, buttonCol1, 7, buttonWidth, 1, handle_save);

        buttons[4].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        buttonCol1 + buttonWidth / 2,
                        3,
                        buttonWidth / 2,
                        1,
                        handle_delete_overlay0);
        buttons[5].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        buttonCol1 + buttonWidth / 2,
                        4,
                        buttonWidth / 2,
                        1,
                        handle_delete_overlay1);
        buttons[6].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        buttonCol1 + buttonWidth / 2,
                        5,
                        buttonWidth / 2,
                        1,
                        handle_delete_overlay2);
        buttons[7].init(("\xf0\x9f\x95\xaf\xef\xb8\x8f"),
                        parent_handle,
                        buttonCol1,
                        2,
                        buttonWidth / 2,
                        1,
                        handle_get_dark_spectrum);
        buttons[8].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        buttonCol1 + buttonWidth / 2,
                        2,
                        buttonWidth / 2,
                        1,
                        handle_delete_dark_spectrum);

        buttons[9].init(("Ref. A"), parent_handle, 0, 11, smallButton, 1, handle_reference_a);
        buttons[10]
            .init(("Ref. B"), parent_handle, smallButton, 11, smallButton, 1, handle_reference_b);
        buttons[11]
            .init(("Reset"), parent_handle, smallButton * 2, 11, smallButton, 1, handle_reset_phase);
        buttons[12]
            .init(("Save"), parent_handle, smallButton * 3, 11, smallButton, 1, handle_save_phase);

        // RGB active
        text_boxes[1].init(parent_handle, textCol1, 2, textWidth, 1);
        text_boxes[2].init(parent_handle, textCol2, 2, textWidth, 1);
        text_boxes[3].init(parent_handle, textCol3, 2, textWidth, 1);

        // RGB overlay0
        text_boxes[4].init(parent_handle, textCol1, 3, textWidth, 1);
        text_boxes[5].init(parent_handle, textCol2, 3, textWidth, 1);
        text_boxes[6].init(parent_handle, textCol3, 3, textWidth, 1);

        // RGB overlay1
        text_boxes[7].init(parent_handle, textCol1, 4, textWidth, 1);
        text_boxes[8].init(parent_handle, textCol2, 4, textWidth, 1);
        text_boxes[9].init(parent_handle, textCol3, 4, textWidth, 1);

        // RGB overlay2
        text_boxes[10].init(parent_handle, textCol1, 5, textWidth, 1);
        text_boxes[11].init(parent_handle, textCol2, 5, textWidth, 1);
        text_boxes[12].init(parent_handle, textCol3, 5, textWidth, 1);

        text_boxes[0].init(parent_handle, textCol3, 7, textWidth, 1);
        text_boxes[0].set_label(-3 * textWidth, 0, "Exposure (ms)");
        text_boxes[0].overwrite_print(std::string("40"));

        text_boxes[13].init(parent_handle, textCol3, 8, textWidth, 1);
        text_boxes[13].set_label(-3 * textWidth, 0, "Pause, count (s, #)");
        text_boxes[14].init(parent_handle, textCol4, 8, textWidth, 1);
        text_boxes[13].overwrite_print(std::string("0"));
        text_boxes[14].overwrite_print(std::string("10"));

        text_boxes[15].init(parent_handle, textCol3, 9, 2, 1);
        text_boxes[15].set_max_characters(6);
        text_boxes[16].init(parent_handle, textCol4, 9, 2, 1);
        text_boxes[16].set_max_characters(6);
        text_boxes[17].init(parent_handle, textCol5, 9, 2, 1);
        text_boxes[17].set_max_characters(6);
        text_boxes[15].set_label(-3 * textWidth, 0, "Freqs. (#, min,max)");
        text_boxes[15].overwrite_print(std::string("2048"));

        text_boxes[18].init(parent_handle, textCol3, 10, 2, 1);
        text_boxes[18].set_max_characters(6);
        text_boxes[19].init(parent_handle, textCol4, 10, 2, 1);
        text_boxes[19].set_max_characters(6);
        text_boxes[20].init(parent_handle, textCol5, 10, 2, 1);
        text_boxes[20].set_max_characters(6);
        text_boxes[18].set_label(-3 * textWidth, 0, "Filter (t0, sig., ord.)");
        text_boxes[18].overwrite_print(std::string("300"));
        text_boxes[19].overwrite_print(std::string("80"));
        text_boxes[20].overwrite_print(std::string("12"));

        file_paths[0].init(parent_handle, 0, 6, 10, 1);
        file_paths[0].set_max_characters(pathChars);
        file_paths[0].overwrite_print(std::format("DefaultOutput.txt"));
        checkboxes[2].init("\xe2\x8c\x9a", parent_handle, buttonCol1 + buttonWidth / 2, 8, 2, 1);
        buttons[17].init(("..."),
                         parent_handle,
                         buttonCol1 + buttonWidth / 2,
                         6,
                         buttonWidth / 2,
                         1,
                         save_path_callback,
                         0);
        text_boxes[22].init(window.parent_handle(4), 3, 0, 2, 1);
        text_boxes[23].init(window.parent_handle(4), 5, 0, 2, 1);
        text_boxes[24].init(window.parent_handle(4), 9, 0, 2, 1);
        text_boxes[25].init(window.parent_handle(4), 11, 0, 2, 1);
        draw_boxes[0].init(window.parent_handle(2), 0, 0, plotWidth, plotHeight);
        draw_boxes[0].set_drawing_function(draw_spectrum);
        checkboxes[1].init(("Log"), window.parent_handle(4), 13, 0, 1, 1);
        checkboxes[3].init(("Average phase"), window.parent_handle(4), 15, 0, 1, 1);
        buttons[13].init(("xlim"), window.parent_handle(4), 2, 0, 1, 1, handle_refresh_request);
        buttons[13].set_tooltip("Apply the entered x limits to the plot. The two text boxes are for "
                               "the upper and lower limits applied to the frequency axis. If they "
                               "are empty, the range will include the whole grid.");
        buttons[14].init(("ylim"), window.parent_handle(4), 7, 0, 1, 1, handle_refresh_request);
        buttons[14].set_tooltip("Apply the entered y limits to the plot. The two text boxes are for "
                               "the upper and lower limits applied to the frequency axis. If they "
                               "are empty, the range will include the whole grid.");
        buttons[15].init(("SVG"), window.parent_handle(4), 1, 0, 1, 1, svg_callback);
        buttons[15].set_tooltip("Generate SVG files of the four line plots, with filenames based on "
                               "the base path set above");
        buttons[16].init("\xe2\x86\x94\xef\xb8\x8f",
                         window.parent_handle(4),
                         0,
                         0,
                         1,
                         1,
                         handle_collapse_panel);
        buttons[16].set_tooltip("Collapse/expand the data entry panel");
        for(int i = 0; i < 18; i++) {
            text_boxes[i].set_max_characters(6);
        }
        pulldowns[1].add_element("Spectrometer live (nm)");
        pulldowns[1].add_element("Spectrometer live (THz)");
        pulldowns[1].add_element("All spectrometers");
        pulldowns[1].add_element("Interferometry: Spectrum");
        pulldowns[1].add_element("Interferometry: Time");
        pulldowns[1].add_element("Interferometry: phase");
        pulldowns[1].add_element("Interferometry: Group delay");
        pulldowns[1].init(parent_handle, 0, 1, 8, 1);
        g_signal_connect(pulldowns[1].element_handle,
                         "notify::selected",
                         G_CALLBACK(drop_down_change_callback),
                         NULL);
        gtk_widget_set_visible(buttons[9].element_handle, false);
        gtk_widget_set_visible(buttons[10].element_handle, false);
        gtk_widget_set_visible(buttons[11].element_handle, false);
        gtk_widget_set_visible(buttons[12].element_handle, false);
        gtk_widget_set_visible(text_boxes[18].label, false);
        gtk_widget_set_visible(text_boxes[18].element_handle, false);
        gtk_widget_set_visible(text_boxes[19].element_handle, false);
        gtk_widget_set_visible(text_boxes[20].element_handle, false);
        console.init(window.parent_handle(1), 0, 0, 1, 1);
        console.c_print("Attached spectrometers:\n");
        initializeSpectrometers();
        pulldowns[0].init(parent_handle, 0, 0, 12, 1);

        for(auto &box : text_boxes) {
            box.vertical_thick();
        }
        for(auto &button : buttons) {
            button.vertical_thick();
        }
        for(auto &p : file_paths) {
            p.vertical_thick();
        }
        for(auto &p : pulldowns) {
            p.vertical_thick();
        }
        window.present();

        draw_boxes[0].queue_draw();
        timeoutID = g_timeout_add(20, G_SOURCE_FUNC(update_display), NULL);
    }

    void initializeSpectrometers() {
        spectrometer_set = std::vector<std::unique_ptr<Spectrometer>>();
        auto ocean_string = open_ocean_spectrometers(spectrometer_set);
        console.c_print("{}", ocean_string);
        for(size_t i = 0; i < spectrometer_set.size(); i++) {
            console.c_print("Device {}: {}\n    serial number: {}\n",
                           i,
                           spectrometer_set[i]->name,
                           spectrometer_set[i]->serial_number);
            std::string new_element = std::format("{}: {}, serial: {}",
                                                  i,
                                                  spectrometer_set[i]->name,
                                                  spectrometer_set[i]->serial_number);
            pulldowns[0].add_element(new_element.c_str());
        }
    }
};
MainGui theGui;

void handle_run_button() {
    if(theGui.running_live()) {
        theGui.stop_live();
    } else {
        theGui.request_live();
    }
}

void svg_callback() {
    theGui.saveSVG = 1;
    theGui.requestPlotUpdate();
}

void handle_refresh_request() { theGui.requestPlotUpdate(); }

void handle_get_overlay0() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).acquire_overlay(0);
}

void drop_down_change_callback() {
    bool visibility = theGui.pulldowns[1].get_value() > 2;
    gtk_widget_set_visible(theGui.buttons[9].element_handle, visibility);
    gtk_widget_set_visible(theGui.buttons[10].element_handle, visibility);
    gtk_widget_set_visible(theGui.buttons[11].element_handle, visibility);
    gtk_widget_set_visible(theGui.buttons[12].element_handle, visibility);
    gtk_widget_set_visible(theGui.text_boxes[18].label, visibility);
    gtk_widget_set_visible(theGui.text_boxes[18].element_handle, visibility);
    gtk_widget_set_visible(theGui.text_boxes[19].element_handle, visibility);
    gtk_widget_set_visible(theGui.text_boxes[20].element_handle, visibility);
}

void handle_get_overlay1() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).acquire_overlay(1);
}

void handle_get_overlay2() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).acquire_overlay(2);
}

void handle_collapse_panel() { theGui.window.toggle_settings_panel(); }

void handle_delete_overlay0() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).delete_overlay(0);
}

void handle_delete_overlay1() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).delete_overlay(1);
}

void handle_delete_overlay2() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).delete_overlay(2);
}

void handle_get_dark_spectrum() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).acquire_dark_spectrum();
}

void handle_delete_dark_spectrum() {
    (*spectrometer_set[theGui.pulldowns[0].get_value()]).disable_dark_spectrum();
}

void acquisition_thread(int active_spectrometer,
                        size_t N,
                        double integration_time,
                        double seconds_to_wait,
                        std::string path,
                        bool timestamp) {
    theBatch.acquire_batch(N,
                           integration_time,
                           seconds_to_wait,
                           (*spectrometer_set[active_spectrometer]));
    theBatch.save(path, timestamp);
    theGui.console.t_print("Finished writing {}!\n", path);
}

void handle_save() {
    theGui.stop_live();
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if((*spectrometer_set[active_spectrometer]).check_lock())
        return;
    std::string path;
    theGui.file_paths[0].copy_buffer(path);
    bool timestamp = theGui.checkboxes[2].is_checked();
    size_t N = (size_t)theGui.text_boxes[14].value_double();
    double integration_time = theGui.text_boxes[0].value_double();
    double waitTime = theGui.text_boxes[13].value_double();
    std::thread(acquisition_thread,
                active_spectrometer,
                N,
                integration_time,
                waitTime,
                path,
                timestamp)
        .detach();
}

void handle_save_phase() {
    std::string path;
    theGui.file_paths[0].copy_buffer(path);
    bool timestamp = theGui.checkboxes[2].is_checked();
    the_interference_controller.save(path, timestamp);
}

void reference_a_acquisition_thread(int active_spectrometer,
                                    size_t N,
                                    double integration_time,
                                    double seconds_to_wait) {
    the_interference_controller.acquire_reference_a(theBatch,
                                                    N,
                                                    integration_time,
                                                    seconds_to_wait,
                                                    (*spectrometer_set[active_spectrometer]));
}

void reference_b_acquisition_thread(int active_spectrometer,
                                    size_t N,
                                    double integration_time,
                                    double seconds_to_wait) {
    the_interference_controller.acquire_reference_b(theBatch,
                                                  N,
                                                  integration_time,
                                                  seconds_to_wait,
                                                  (*spectrometer_set[active_spectrometer]));
}

void handle_reference_a() {
    handle_reset_controller();
    if(!the_interference_controller.check_configuration_status()) {
        theGui.console.c_print("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stop_live();
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if((*spectrometer_set[active_spectrometer]).check_lock())
        return;
    size_t N = (size_t)theGui.text_boxes[14].value_double();
    double integration_time = theGui.text_boxes[0].value_double();
    double waitTime = theGui.text_boxes[13].value_double();
    std::thread(reference_a_acquisition_thread, active_spectrometer, N, integration_time, waitTime)
        .detach();
}

void handle_reference_b() {
    handle_reset_controller();
    if(!the_interference_controller.check_configuration_status()) {
        theGui.console.c_print("Reset controller with frequencies first.\n");
        return;
    }
    theGui.stop_live();
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if((*spectrometer_set[active_spectrometer]).check_lock())
        return;
    size_t N = (size_t)theGui.text_boxes[14].value_double();
    double integration_time = theGui.text_boxes[0].value_double();
    double waitTime = theGui.text_boxes[13].value_double();
    std::thread(reference_b_acquisition_thread, active_spectrometer, N, integration_time, waitTime)
        .detach();
}

void handle_reset_controller() {
    size_t num_freqs = static_cast<size_t>(theGui.text_boxes[15].value_double());
    if(num_freqs < 2)
        return;
    double fMin = theGui.text_boxes[16].value_double();
    double fMax = theGui.text_boxes[17].value_double();
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(fMax == 0.0) {
        fMax = (*spectrometer_set[active_spectrometer]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if(fMin == 0.0) {
        fMin = (*spectrometer_set[active_spectrometer])
                   .wavelengths()[(*spectrometer_set[active_spectrometer]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    the_interference_controller.reset_frequencies(num_freqs, fMin, fMax);
}

void handle_reset_phase() { the_interference_controller.reset_phase(); }

void draw_spectrum(GtkDrawingArea *area, cairo_t *cr, int width, int height, gpointer data) {
    switch(theGui.pulldowns[1].get_value()) {
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

    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    if(active_spectrometer != -1 &&
       (active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).initialized()) {
        return;
    }

    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceYmin = true;
        forceYmax = yMax != 0.0;
    }
    if(logPlot && !forceYmin) {
        forceYmin = true;
        yMin = 1e-1;
    }

    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_Spectrum.svg");
        sPlot.SVGPath = svgPath;
    }

    if(theGui.running_live() && (active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).check_lock()) {
        (*spectrometer_set[active_spectrometer])
            .set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].value_double()));
        (*spectrometer_set[active_spectrometer]).acquire_single();
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataX = (*spectrometer_set[active_spectrometer]).wavelengths();
    sPlot.hasDataX = true;
    sPlot.data = (*spectrometer_set[active_spectrometer]).data();
    sPlot.Npts = (*spectrometer_set[active_spectrometer]).size();
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceYmin;
    sPlot.forceYmax = forceYmax;
    if(forceYmax)
        sPlot.forcedYmax = yMax;
    if(forceYmin)
        sPlot.forcedYmin = yMin;
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
    if((*spectrometer_set[active_spectrometer]).get_overlay_count() > 0) {
        int firstAdded = -1;
        sPlot.ExtraLines = (*spectrometer_set[active_spectrometer]).get_overlay_count();
        if((*spectrometer_set[active_spectrometer]).has_overlay_0) {
            sPlot.data2 = (*spectrometer_set[active_spectrometer]).get_overlay(0);
            firstAdded = 0;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_1) {
            sPlot.data2 = (*spectrometer_set[active_spectrometer]).get_overlay(1);
            firstAdded = 1;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_2) {
            sPlot.data2 = (*spectrometer_set[active_spectrometer]).get_overlay(2);
            firstAdded = 2;
        }
        if(sPlot.ExtraLines > 1 && firstAdded != 1 &&
           (*spectrometer_set[active_spectrometer]).has_overlay_1)
            sPlot.data3 = (*spectrometer_set[active_spectrometer]).get_overlay(1);
        else if(sPlot.ExtraLines > 1 && firstAdded != 2 &&
                (*spectrometer_set[active_spectrometer]).has_overlay_2)
            sPlot.data3 = (*spectrometer_set[active_spectrometer]).get_overlay(2);

        if(sPlot.ExtraLines > 2)
            sPlot.data4 = (*spectrometer_set[active_spectrometer]).get_overlay(2);
    }
    sPlot.plot(cr);
}

void draw_spectrum_frequency(GtkDrawingArea *area,
                             cairo_t *cr,
                             int width,
                             int height,
                             gpointer data) {
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    if((active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).initialized()) {
        return;
    }

    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceYmin = true;
        forceYmax = yMax != 0.0;
    }
    if(logPlot && !forceYmin) {
        forceYmin = true;
        yMin = 1e-1;
    }

    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t num_freqs = static_cast<size_t>(theGui.text_boxes[15].value_double());
    if(num_freqs < 2)
        return;
    double fMin = theGui.text_boxes[16].value_double();
    double fMax = theGui.text_boxes[17].value_double();

    if(fMax == 0.0) {
        fMax = (*spectrometer_set[active_spectrometer]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
    }
    if(fMin == 0.0) {
        fMin = (*spectrometer_set[active_spectrometer])
                   .wavelengths()[(*spectrometer_set[active_spectrometer]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
    }
    double dF = (fMax - fMin) / static_cast<double>(num_freqs - 1);
    std::vector<double> frequencies(num_freqs);
    for(size_t i = 0; i < num_freqs; i++) {
        frequencies[i] = fMin + static_cast<double>(i) * dF;
    }
    std::vector<double> liveSpectrum;
    if(theGui.running_live() && (active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).check_lock() && num_freqs > 0) {
        (*spectrometer_set[active_spectrometer])
            .set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].value_double()));
        liveSpectrum =
            (*spectrometer_set[active_spectrometer]).acquire_single_frequency(frequencies);
    } else
        return;

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
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
    if(forceYmax)
        sPlot.forcedYmax = yMax;
    if(forceYmin)
        sPlot.forcedYmin = yMin;
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
    if((*spectrometer_set[active_spectrometer]).get_overlay_count() > 0) {
        int firstAdded = -1;
        sPlot.ExtraLines = (*spectrometer_set[active_spectrometer]).get_overlay_count();
        if((*spectrometer_set[active_spectrometer]).has_overlay_0) {
            sPlot.data2 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(0, frequencies);
            firstAdded = 0;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_1) {
            sPlot.data2 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(1, frequencies);
            firstAdded = 1;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_2) {
            sPlot.data2 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(2, frequencies);
            firstAdded = 2;
        }
        if(sPlot.ExtraLines > 1 && firstAdded != 1 &&
           (*spectrometer_set[active_spectrometer]).has_overlay_1)
            sPlot.data3 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(1, frequencies);
        else if(sPlot.ExtraLines > 1 && firstAdded != 2 &&
                (*spectrometer_set[active_spectrometer]).has_overlay_2)
            sPlot.data3 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(2, frequencies);

        if(sPlot.ExtraLines > 2)
            sPlot.data4 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(2, frequencies);
    }
    sPlot.plot(cr);
}

void draw_spectra_frequency(GtkDrawingArea *area,
                            cairo_t *cr,
                            int width,
                            int height,
                            gpointer data) {
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceYmin = false;
    bool forceYmax = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceYmin = true;
        forceYmax = yMax != 0.0;
    }
    if(logPlot && !forceYmin) {
        forceYmin = true;
        yMin = 1e-1;
    }
    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_SpectrumTHz.svg");
        sPlot.SVGPath = svgPath;
    }
    size_t num_freqs = static_cast<size_t>(theGui.text_boxes[15].value_double());
    if(num_freqs < 2)
        return;
    double fMin = theGui.text_boxes[16].value_double();
    double fMax = theGui.text_boxes[17].value_double();
    if(fMax == 0.0) {
        fMax = (*spectrometer_set[0]).wavelengths()[0];
        fMax = 1e-3 * lightC<double>() / fMax;
        for(size_t i = 1; i < spectrometer_set.size(); i++) {
            fMax = maxN(fMax, 1e-3 * lightC<double>() / (*spectrometer_set[i]).wavelengths()[0]);
        }
    }
    if(fMin == 0.0) {
        fMin = (*spectrometer_set[0]).wavelengths()[(*spectrometer_set[0]).size() - 1];
        fMin = 1e-3 * lightC<double>() / fMin;
        for(size_t i = 1; i < spectrometer_set.size(); i++) {
            fMin =
                minN(fMin,
                     1e-3 * lightC<double>() /
                         (*spectrometer_set[i]).wavelengths()[(*spectrometer_set[0]).size() - 1]);
        }
    }

    double dF = (fMax - fMin) / static_cast<double>(num_freqs - 1);
    std::vector<double> frequencies(num_freqs);
    for(size_t i = 0; i < num_freqs; i++) {
        frequencies[i] = fMin + static_cast<double>(i) * dF;
    }
    std::vector<std::vector<double>> liveSpectra(spectrometer_set.size());
    if(theGui.running_live() && num_freqs > 0) {
        for(size_t i = 0; i < spectrometer_set.size(); i++) {
            liveSpectra[i] = (*spectrometer_set[i]).acquire_single_frequency(frequencies);
        }
    } else
        return;

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
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
    if(forceYmax)
        sPlot.forcedYmax = yMax;
    if(forceYmin)
        sPlot.forcedYmin = yMin;
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
    sPlot.ExtraLines = static_cast<int>(liveSpectra.size()) - 1;
    if(liveSpectra.size() > 1)
        sPlot.data2 = liveSpectra[1].data();
    if(liveSpectra.size() > 2)
        sPlot.data3 = liveSpectra[2].data();
    if(liveSpectra.size() > 3)
        sPlot.data4 = liveSpectra[3].data();
    sPlot.plot(cr);
}

void draw_interference_spectrum(GtkDrawingArea *area,
                                cairo_t *cr,
                                int width,
                                int height,
                                gpointer data) {
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    if(!the_interference_controller.check_configuration_status())
        handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_InterferenceF.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
    }

    if(theGui.running_live() && !(*spectrometer_set[theGui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[theGui.pulldowns[0].get_value()])
            .set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_interferogram(
            (*spectrometer_set[theGui.pulldowns[0].get_value()]));
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
    if(forceY)
        sPlot.forcedYmin = yMin;
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
    if(the_interference_controller.get_reference_a_vector().size() ==
       the_interference_controller.get_num_freqs()) {
        sPlot.data2 = the_interference_controller.get_reference_a();
        sPlot.ExtraLines = 1;
        if(the_interference_controller.get_reference_b_vector().size() ==
           the_interference_controller.get_num_freqs()) {
            sPlot.data3 = the_interference_controller.get_reference_b();
            sPlot.ExtraLines = 2;
        }
    }
    sPlot.plot(cr);
}

void draw_interference_spectrum_time(GtkDrawingArea *area,
                                     cairo_t *cr,
                                     int width,
                                     int height,
                                     gpointer data) {
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_InterferenceF.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
    }
    double t0 = theGui.text_boxes[18].value_double();
    double sigma = theGui.text_boxes[19].value_double();
    double ord = theGui.text_boxes[20].value_double();
    if(theGui.running_live() && !(*spectrometer_set[theGui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[theGui.pulldowns[0].get_value()])
            .set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_interferogram(
            (*spectrometer_set[theGui.pulldowns[0].get_value()]));
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
    if(forceY)
        sPlot.forcedYmin = yMin;
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
    sPlot.data2 = the_interference_controller.get_time_reference_a();
    sPlot.data3 = the_interference_controller.get_time_reference_b();
    sPlot.data4 = the_interference_controller.get_time_filter();
    sPlot.ExtraLines = 3;
    sPlot.plot(cr);
}

void draw_interference_phase(GtkDrawingArea *area,
                             cairo_t *cr,
                             int width,
                             int height,
                             gpointer data) {
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_InterferencePhase.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
    }

    if(theGui.running_live() && !(*spectrometer_set[theGui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[theGui.pulldowns[0].get_value()])
            .set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_phase(
            (*spectrometer_set[theGui.pulldowns[0].get_value()]));
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
    if(forceY)
        sPlot.forcedYmin = yMin;
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

void draw_interference_group_delay(GtkDrawingArea *area,
                                   cairo_t *cr,
                                   int width,
                                   int height,
                                   gpointer data) {
    int active_spectrometer = theGui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if(saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if(theGui.checkboxes[1].is_checked()) {
        logPlot = true;
    }

    bool forceX = false;
    double xMin = theGui.text_boxes[22].value_double();
    double xMax = theGui.text_boxes[23].value_double();
    if(xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.text_boxes[24].value_double();
    double yMax = theGui.text_boxes[25].value_double();
    if(yMin != yMax && yMax > yMin) {
        forceY = true;
    }

    if(saveSVG) {
        std::string svgPath;
        theGui.file_paths[0].copy_buffer(svgPath);
        svgPath.append("_InterferencePhase.svg");
        sPlot.SVGPath = svgPath;
    }

    LweColor mainColor(0.5, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[1].value_double() + theGui.text_boxes[2].value_double() +
               theGui.text_boxes[3].value_double())) {
        mainColor = LweColor(theGui.text_boxes[1].value_double(),
                             theGui.text_boxes[2].value_double(),
                             theGui.text_boxes[3].value_double(),
                             1);
    }

    LweColor overLay0Color(1, 0, 1, 1);
    if(0.0 != (theGui.text_boxes[4].value_double() + theGui.text_boxes[5].value_double() +
               theGui.text_boxes[6].value_double())) {
        overLay0Color = LweColor(theGui.text_boxes[4].value_double(),
                                 theGui.text_boxes[5].value_double(),
                                 theGui.text_boxes[6].value_double(),
                                 1);
    }

    LweColor overLay1Color(1, 0.5, 0, 1);
    if(0.0 != (theGui.text_boxes[7].value_double() + theGui.text_boxes[8].value_double() +
               theGui.text_boxes[9].value_double())) {
        overLay1Color = LweColor(theGui.text_boxes[7].value_double(),
                                 theGui.text_boxes[8].value_double(),
                                 theGui.text_boxes[9].value_double(),
                                 1);
    }

    LweColor overLay2Color(0, 1, 1, 1);
    if(0.0 != (theGui.text_boxes[10].value_double() + theGui.text_boxes[11].value_double() +
               theGui.text_boxes[12].value_double())) {
        overLay2Color = LweColor(theGui.text_boxes[10].value_double(),
                                 theGui.text_boxes[11].value_double(),
                                 theGui.text_boxes[12].value_double(),
                                 1);
    }

    if(theGui.running_live() && !(*spectrometer_set[theGui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[theGui.pulldowns[0].get_value()])
            .set_integration_time((unsigned long)round(1000 * theGui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(theGui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_phase(
            (*spectrometer_set[theGui.pulldowns[0].get_value()]));
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
    if(forceY)
        sPlot.forcedYmin = yMin;
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
    theGui.console.update_from_buffer();
    theGui.apply_update();
    return true;
}

void save_path_callback() {
#ifdef __APPLE__
    theGui.path_buffer = pathFromAppleSaveDialog();
#endif
    theGui.requestSavePathUpdate();
}

static void activate(GtkApplication *app, gpointer user_data) {
#if defined __linux__ || defined __APPLE__
    setlocale(LC_NUMERIC, "en_US.UTF-8");
#else
    setlocale(LC_NUMERIC, "en_US");
#endif
    theGui.activate(app);
}

int main(int argc, char **argv) {
    GtkApplication *app =
        gtk_application_new("io.github.NickKarpowicz.Scarab", (GApplicationFlags)0);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    int status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    return status;
}

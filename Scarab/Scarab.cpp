#include "Scarab.hpp"
#include "batch_acquisition.hpp"
#include "gtk_wrappers.hpp"
#include "interferometry.hpp"
#include "ocean_spectrometer.hpp"
#include "spectrometer.hpp"
#include <thread>

std::vector<std::unique_ptr<Spectrometer>> spectrometer_set;
BatchAcquisition the_batch;
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
    int save_svg = 0;
    bool loaded_defaults = false;
    unsigned int timeout_id = 0;
    void request_plot_update() { queue_update = true; }

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

    void request_interface_values_update() { queue_interface_values_update = true; }
    void request_save_path_update() {
#ifndef __APPLE__
        path_buffer = std::string("?LWE_LOADING??");
        path_from_save_dialog(path_buffer, "txt", "Text (.txt)");
#endif
        queue_save_path_update = true;
    }
    void activate(GtkApplication *app) {
        int button_width = 4;
        int small_button = 3;
        int text_width = 2;
        int plot_width = 12;
        int plot_height = 6;
        int path_chars = 42;
        int text_col1 = text_width;
        int text_col2 = 2 * text_width;
        int text_col3 = 3 * text_width;
        int text_col4 = 4 * text_width;
        int text_col5 = 5 * text_width;
        int button_col1 = 8;

        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", true, NULL);
        window.init(app, "Spectrometer Control and Recording for Attosecond Beamlines", 1400, 800);
        GtkWidget *parent_handle = window.parent_handle();

        buttons[0].init(("Run"), parent_handle, button_col1, 1, button_width, 1, handle_run_button);
        buttons[1].init(("\xf0\x9f\x93\x88"),
                        parent_handle,
                        button_col1,
                        3,
                        button_width / 2,
                        1,
                        handle_get_overlay0);
        buttons[2].init(("\xf0\x9f\x93\x88"),
                        parent_handle,
                        button_col1,
                        4,
                        button_width / 2,
                        1,
                        handle_get_overlay1);
        buttons[3].init(("\xf0\x9f\x93\x88"),
                        parent_handle,
                        button_col1,
                        5,
                        button_width / 2,
                        1,
                        handle_get_overlay2);
        buttons[18].init(("Acquire"), parent_handle, button_col1, 7, button_width, 1, handle_save);

        buttons[4].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        button_col1 + button_width / 2,
                        3,
                        button_width / 2,
                        1,
                        handle_delete_overlay0);
        buttons[5].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        button_col1 + button_width / 2,
                        4,
                        button_width / 2,
                        1,
                        handle_delete_overlay1);
        buttons[6].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        button_col1 + button_width / 2,
                        5,
                        button_width / 2,
                        1,
                        handle_delete_overlay2);
        buttons[7].init(("\xf0\x9f\x95\xaf\xef\xb8\x8f"),
                        parent_handle,
                        button_col1,
                        2,
                        button_width / 2,
                        1,
                        handle_get_dark_spectrum);
        buttons[8].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"),
                        parent_handle,
                        button_col1 + button_width / 2,
                        2,
                        button_width / 2,
                        1,
                        handle_delete_dark_spectrum);

        buttons[9].init(("Ref. A"), parent_handle, 0, 11, small_button, 1, handle_reference_a);
        buttons[10]
            .init(("Ref. B"), parent_handle, small_button, 11, small_button, 1, handle_reference_b);
        buttons[11].init(("Reset"),
                         parent_handle,
                         small_button * 2,
                         11,
                         small_button,
                         1,
                         handle_reset_phase);
        buttons[12].init(("Save"),
                         parent_handle,
                         small_button * 3,
                         11,
                         small_button,
                         1,
                         handle_save_phase);

        // RGB active
        text_boxes[1].init(parent_handle, text_col1, 2, text_width, 1);
        text_boxes[2].init(parent_handle, text_col2, 2, text_width, 1);
        text_boxes[3].init(parent_handle, text_col3, 2, text_width, 1);

        // RGB overlay0
        text_boxes[4].init(parent_handle, text_col1, 3, text_width, 1);
        text_boxes[5].init(parent_handle, text_col2, 3, text_width, 1);
        text_boxes[6].init(parent_handle, text_col3, 3, text_width, 1);

        // RGB overlay1
        text_boxes[7].init(parent_handle, text_col1, 4, text_width, 1);
        text_boxes[8].init(parent_handle, text_col2, 4, text_width, 1);
        text_boxes[9].init(parent_handle, text_col3, 4, text_width, 1);

        // RGB overlay2
        text_boxes[10].init(parent_handle, text_col1, 5, text_width, 1);
        text_boxes[11].init(parent_handle, text_col2, 5, text_width, 1);
        text_boxes[12].init(parent_handle, text_col3, 5, text_width, 1);

        text_boxes[0].init(parent_handle, text_col3, 7, text_width, 1);
        text_boxes[0].set_label(-3 * text_width, 0, "Exposure (ms)");
        text_boxes[0].overwrite_print(std::string("40"));

        text_boxes[13].init(parent_handle, text_col3, 8, text_width, 1);
        text_boxes[13].set_label(-3 * text_width, 0, "Pause, count (s, #)");
        text_boxes[14].init(parent_handle, text_col4, 8, text_width, 1);
        text_boxes[13].overwrite_print(std::string("0"));
        text_boxes[14].overwrite_print(std::string("10"));

        text_boxes[15].init(parent_handle, text_col3, 9, 2, 1);
        text_boxes[15].set_max_characters(6);
        text_boxes[16].init(parent_handle, text_col4, 9, 2, 1);
        text_boxes[16].set_max_characters(6);
        text_boxes[17].init(parent_handle, text_col5, 9, 2, 1);
        text_boxes[17].set_max_characters(6);
        text_boxes[15].set_label(-3 * text_width, 0, "Freqs. (#, min,max)");
        text_boxes[15].overwrite_print(std::string("2048"));

        text_boxes[18].init(parent_handle, text_col3, 10, 2, 1);
        text_boxes[18].set_max_characters(6);
        text_boxes[19].init(parent_handle, text_col4, 10, 2, 1);
        text_boxes[19].set_max_characters(6);
        text_boxes[20].init(parent_handle, text_col5, 10, 2, 1);
        text_boxes[20].set_max_characters(6);
        text_boxes[18].set_label(-3 * text_width, 0, "Filter (t0, sig., ord.)");
        text_boxes[18].overwrite_print(std::string("300"));
        text_boxes[19].overwrite_print(std::string("80"));
        text_boxes[20].overwrite_print(std::string("12"));

        file_paths[0].init(parent_handle, 0, 6, 10, 1);
        file_paths[0].set_max_characters(path_chars);
        file_paths[0].overwrite_print(std::format("DefaultOutput.txt"));
        checkboxes[2].init("\xe2\x8c\x9a", parent_handle, button_col1 + button_width / 2, 8, 2, 1);
        buttons[17].init(("..."),
                         parent_handle,
                         button_col1 + button_width / 2,
                         6,
                         button_width / 2,
                         1,
                         save_path_callback,
                         0);
        text_boxes[22].init(window.parent_handle(4), 3, 0, 2, 1);
        text_boxes[23].init(window.parent_handle(4), 5, 0, 2, 1);
        text_boxes[24].init(window.parent_handle(4), 9, 0, 2, 1);
        text_boxes[25].init(window.parent_handle(4), 11, 0, 2, 1);
        draw_boxes[0].init(window.parent_handle(2), 0, 0, plot_width, plot_height);
        draw_boxes[0].set_drawing_function(draw_spectrum);
        checkboxes[1].init(("Log"), window.parent_handle(4), 13, 0, 1, 1);
        checkboxes[3].init(("Average phase"), window.parent_handle(4), 15, 0, 1, 1);
        buttons[13].init(("xlim"), window.parent_handle(4), 2, 0, 1, 1, handle_refresh_request);
        buttons[13].set_tooltip(
            "Apply the entered x limits to the plot. The two text boxes are for "
            "the upper and lower limits applied to the frequency axis. If they "
            "are empty, the range will include the whole grid.");
        buttons[14].init(("ylim"), window.parent_handle(4), 7, 0, 1, 1, handle_refresh_request);
        buttons[14].set_tooltip(
            "Apply the entered y limits to the plot. The two text boxes are for "
            "the upper and lower limits applied to the frequency axis. If they "
            "are empty, the range will include the whole grid.");
        buttons[15].init(("SVG"), window.parent_handle(4), 1, 0, 1, 1, svg_callback);
        buttons[15].set_tooltip(
            "Generate SVG files of the four line plots, with filenames based on "
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
        initialize_spectrometers();
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
        timeout_id = g_timeout_add(20, G_SOURCE_FUNC(update_display), NULL);
    }

    void initialize_spectrometers() {
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
MainGui the_gui;

void handle_run_button() {
    if(the_gui.running_live()) {
        the_gui.stop_live();
    } else {
        the_gui.request_live();
    }
}

void svg_callback() {
    the_gui.save_svg = 1;
    the_gui.request_plot_update();
}

void handle_refresh_request() { the_gui.request_plot_update(); }

void handle_get_overlay0() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).acquire_overlay(0);
}

void drop_down_change_callback() {
    bool visibility = the_gui.pulldowns[1].get_value() > 2;
    gtk_widget_set_visible(the_gui.buttons[9].element_handle, visibility);
    gtk_widget_set_visible(the_gui.buttons[10].element_handle, visibility);
    gtk_widget_set_visible(the_gui.buttons[11].element_handle, visibility);
    gtk_widget_set_visible(the_gui.buttons[12].element_handle, visibility);
    gtk_widget_set_visible(the_gui.text_boxes[18].label, visibility);
    gtk_widget_set_visible(the_gui.text_boxes[18].element_handle, visibility);
    gtk_widget_set_visible(the_gui.text_boxes[19].element_handle, visibility);
    gtk_widget_set_visible(the_gui.text_boxes[20].element_handle, visibility);
}

void handle_get_overlay1() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).acquire_overlay(1);
}

void handle_get_overlay2() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).acquire_overlay(2);
}

void handle_collapse_panel() { the_gui.window.toggle_settings_panel(); }

void handle_delete_overlay0() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).delete_overlay(0);
}

void handle_delete_overlay1() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).delete_overlay(1);
}

void handle_delete_overlay2() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).delete_overlay(2);
}

void handle_get_dark_spectrum() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).acquire_dark_spectrum();
}

void handle_delete_dark_spectrum() {
    (*spectrometer_set[the_gui.pulldowns[0].get_value()]).disable_dark_spectrum();
}

void acquisition_thread(int active_spectrometer,
                        size_t n,
                        double integration_time,
                        double seconds_to_wait,
                        std::string path,
                        bool timestamp) {
    the_batch.acquire_batch(n,
                            integration_time,
                            seconds_to_wait,
                            (*spectrometer_set[active_spectrometer]));
    the_batch.save(path, timestamp);
    the_gui.console.t_print("Finished writing {}!\n", path);
}

void handle_save() {
    the_gui.stop_live();
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if((*spectrometer_set[active_spectrometer]).check_lock())
        return;
    std::string path;
    the_gui.file_paths[0].copy_buffer(path);
    bool timestamp = the_gui.checkboxes[2].is_checked();
    size_t n = (size_t)the_gui.text_boxes[14].value_double();
    double integration_time = the_gui.text_boxes[0].value_double();
    double wait_time = the_gui.text_boxes[13].value_double();
    std::thread(acquisition_thread,
                active_spectrometer,
                n,
                integration_time,
                wait_time,
                path,
                timestamp)
        .detach();
}

void handle_save_phase() {
    std::string path;
    the_gui.file_paths[0].copy_buffer(path);
    bool timestamp = the_gui.checkboxes[2].is_checked();
    the_interference_controller.save(path, timestamp);
}

void reference_a_acquisition_thread(int active_spectrometer,
                                    size_t n,
                                    double integration_time,
                                    double seconds_to_wait) {
    the_interference_controller.acquire_reference_a(the_batch,
                                                    n,
                                                    integration_time,
                                                    seconds_to_wait,
                                                    (*spectrometer_set[active_spectrometer]));
}

void reference_b_acquisition_thread(int active_spectrometer,
                                    size_t n,
                                    double integration_time,
                                    double seconds_to_wait) {
    the_interference_controller.acquire_reference_b(the_batch,
                                                    n,
                                                    integration_time,
                                                    seconds_to_wait,
                                                    (*spectrometer_set[active_spectrometer]));
}

void handle_reference_a() {
    handle_reset_controller();
    if(!the_interference_controller.check_configuration_status()) {
        the_gui.console.c_print("Reset controller with frequencies first.\n");
        return;
    }
    the_gui.stop_live();
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if((*spectrometer_set[active_spectrometer]).check_lock())
        return;
    size_t n = (size_t)the_gui.text_boxes[14].value_double();
    double integration_time = the_gui.text_boxes[0].value_double();
    double wait_time = the_gui.text_boxes[13].value_double();
    std::thread(reference_a_acquisition_thread, active_spectrometer, n, integration_time, wait_time)
        .detach();
}

void handle_reference_b() {
    handle_reset_controller();
    if(!the_interference_controller.check_configuration_status()) {
        the_gui.console.c_print("Reset controller with frequencies first.\n");
        return;
    }
    the_gui.stop_live();
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if((*spectrometer_set[active_spectrometer]).check_lock())
        return;
    size_t n = (size_t)the_gui.text_boxes[14].value_double();
    double integration_time = the_gui.text_boxes[0].value_double();
    double wait_time = the_gui.text_boxes[13].value_double();
    std::thread(reference_b_acquisition_thread, active_spectrometer, n, integration_time, wait_time)
        .detach();
}

void handle_reset_controller() {
    size_t num_freqs = static_cast<size_t>(the_gui.text_boxes[15].value_double());
    if(num_freqs < 2)
        return;
    double f_min = the_gui.text_boxes[16].value_double();
    double f_max = the_gui.text_boxes[17].value_double();
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(f_max == 0.0) {
        f_max = (*spectrometer_set[active_spectrometer]).wavelengths()[0];
        f_max = 1e-3 * lightC<double>() / f_max;
    }
    if(f_min == 0.0) {
        f_min = (*spectrometer_set[active_spectrometer])
                    .wavelengths()[(*spectrometer_set[active_spectrometer]).size() - 1];
        f_min = 1e-3 * lightC<double>() / f_min;
    }
    the_interference_controller.reset_frequencies(num_freqs, f_min, f_max);
}

void handle_reset_phase() { the_interference_controller.reset_phase(); }

void draw_spectrum(GtkDrawingArea *area, cairo_t *cr, int width, int height, gpointer data) {
    switch(the_gui.pulldowns[1].get_value()) {
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

    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    if(active_spectrometer != -1 &&
       (active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).initialized()) {
        return;
    }

    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_ymin = false;
    bool force_ymax = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_ymin = true;
        force_ymax = y_max != 0.0;
    }
    if(log_plot && !force_ymin) {
        force_ymin = true;
        y_min = 1e-1;
    }

    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_Spectrum.svg");
        s_plot.SVGPath = svg_path;
    }

    if(the_gui.running_live() &&
       (active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).check_lock()) {
        (*spectrometer_set[active_spectrometer])
            .set_integration_time(
                (unsigned long)round(1000 * the_gui.text_boxes[0].value_double()));
        (*spectrometer_set[active_spectrometer]).acquire_single();
    }

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }

    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = (*spectrometer_set[active_spectrometer]).wavelengths();
    s_plot.hasDataX = true;
    s_plot.data = (*spectrometer_set[active_spectrometer]).data();
    s_plot.Npts = (*spectrometer_set[active_spectrometer]).size();
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_ymin;
    s_plot.forceYmax = force_ymax;
    if(force_ymax)
        s_plot.forcedYmax = y_max;
    if(force_ymin)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Wavelength (nm)";
    s_plot.yLabel = "Spectrum (counts)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    if((*spectrometer_set[active_spectrometer]).get_overlay_count() > 0) {
        int first_added = -1;
        s_plot.ExtraLines = (*spectrometer_set[active_spectrometer]).get_overlay_count();
        if((*spectrometer_set[active_spectrometer]).has_overlay_0) {
            s_plot.data2 = (*spectrometer_set[active_spectrometer]).get_overlay(0);
            first_added = 0;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_1) {
            s_plot.data2 = (*spectrometer_set[active_spectrometer]).get_overlay(1);
            first_added = 1;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_2) {
            s_plot.data2 = (*spectrometer_set[active_spectrometer]).get_overlay(2);
            first_added = 2;
        }
        if(s_plot.ExtraLines > 1 && first_added != 1 &&
           (*spectrometer_set[active_spectrometer]).has_overlay_1)
            s_plot.data3 = (*spectrometer_set[active_spectrometer]).get_overlay(1);
        else if(s_plot.ExtraLines > 1 && first_added != 2 &&
                (*spectrometer_set[active_spectrometer]).has_overlay_2)
            s_plot.data3 = (*spectrometer_set[active_spectrometer]).get_overlay(2);

        if(s_plot.ExtraLines > 2)
            s_plot.data4 = (*spectrometer_set[active_spectrometer]).get_overlay(2);
    }
    s_plot.plot(cr);
}

void draw_spectrum_frequency(GtkDrawingArea *area,
                             cairo_t *cr,
                             int width,
                             int height,
                             gpointer data) {
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    if((active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).initialized()) {
        return;
    }

    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_ymin = false;
    bool force_ymax = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_ymin = true;
        force_ymax = y_max != 0.0;
    }
    if(log_plot && !force_ymin) {
        force_ymin = true;
        y_min = 1e-1;
    }

    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_SpectrumTHz.svg");
        s_plot.SVGPath = svg_path;
    }
    size_t num_freqs = static_cast<size_t>(the_gui.text_boxes[15].value_double());
    if(num_freqs < 2)
        return;
    double f_min = the_gui.text_boxes[16].value_double();
    double f_max = the_gui.text_boxes[17].value_double();

    if(f_max == 0.0) {
        f_max = (*spectrometer_set[active_spectrometer]).wavelengths()[0];
        f_max = 1e-3 * lightC<double>() / f_max;
    }
    if(f_min == 0.0) {
        f_min = (*spectrometer_set[active_spectrometer])
                    .wavelengths()[(*spectrometer_set[active_spectrometer]).size() - 1];
        f_min = 1e-3 * lightC<double>() / f_min;
    }
    double d_f = (f_max - f_min) / static_cast<double>(num_freqs - 1);
    std::vector<double> frequencies(num_freqs);
    for(size_t i = 0; i < num_freqs; i++) {
        frequencies[i] = f_min + static_cast<double>(i) * d_f;
    }
    std::vector<double> live_spectrum;
    if(the_gui.running_live() &&
       (active_spectrometer < static_cast<int>(spectrometer_set.size())) &&
       !(*spectrometer_set[active_spectrometer]).check_lock() && num_freqs > 0) {
        (*spectrometer_set[active_spectrometer])
            .set_integration_time(
                (unsigned long)round(1000 * the_gui.text_boxes[0].value_double()));
        live_spectrum =
            (*spectrometer_set[active_spectrometer]).acquire_single_frequency(frequencies);
    } else
        return;

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }

    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = frequencies.data();
    s_plot.hasDataX = true;
    s_plot.data = live_spectrum.data();
    s_plot.Npts = num_freqs;
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_ymin;
    s_plot.forceYmax = force_ymax;
    if(force_ymax)
        s_plot.forcedYmax = y_max;
    if(force_ymin)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Frequency (THz)";
    s_plot.yLabel = "Spectrum (counts)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    if((*spectrometer_set[active_spectrometer]).get_overlay_count() > 0) {
        int first_added = -1;
        s_plot.ExtraLines = (*spectrometer_set[active_spectrometer]).get_overlay_count();
        if((*spectrometer_set[active_spectrometer]).has_overlay_0) {
            s_plot.data2 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(0, frequencies);
            first_added = 0;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_1) {
            s_plot.data2 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(1, frequencies);
            first_added = 1;
        } else if((*spectrometer_set[active_spectrometer]).has_overlay_2) {
            s_plot.data2 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(2, frequencies);
            first_added = 2;
        }
        if(s_plot.ExtraLines > 1 && first_added != 1 &&
           (*spectrometer_set[active_spectrometer]).has_overlay_1)
            s_plot.data3 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(1, frequencies);
        else if(s_plot.ExtraLines > 1 && first_added != 2 &&
                (*spectrometer_set[active_spectrometer]).has_overlay_2)
            s_plot.data3 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(2, frequencies);

        if(s_plot.ExtraLines > 2)
            s_plot.data4 =
                (*spectrometer_set[active_spectrometer]).get_overlay_frequency(2, frequencies);
    }
    s_plot.plot(cr);
}

void draw_spectra_frequency(GtkDrawingArea *area,
                            cairo_t *cr,
                            int width,
                            int height,
                            gpointer data) {
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_ymin = false;
    bool force_ymax = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_ymin = true;
        force_ymax = y_max != 0.0;
    }
    if(log_plot && !force_ymin) {
        force_ymin = true;
        y_min = 1e-1;
    }
    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_SpectrumTHz.svg");
        s_plot.SVGPath = svg_path;
    }
    size_t num_freqs = static_cast<size_t>(the_gui.text_boxes[15].value_double());
    if(num_freqs < 2)
        return;
    double f_min = the_gui.text_boxes[16].value_double();
    double f_max = the_gui.text_boxes[17].value_double();
    if(f_max == 0.0) {
        f_max = (*spectrometer_set[0]).wavelengths()[0];
        f_max = 1e-3 * lightC<double>() / f_max;
        for(size_t i = 1; i < spectrometer_set.size(); i++) {
            f_max = maxN(f_max, 1e-3 * lightC<double>() / (*spectrometer_set[i]).wavelengths()[0]);
        }
    }
    if(f_min == 0.0) {
        f_min = (*spectrometer_set[0]).wavelengths()[(*spectrometer_set[0]).size() - 1];
        f_min = 1e-3 * lightC<double>() / f_min;
        for(size_t i = 1; i < spectrometer_set.size(); i++) {
            f_min =
                minN(f_min,
                     1e-3 * lightC<double>() /
                         (*spectrometer_set[i]).wavelengths()[(*spectrometer_set[0]).size() - 1]);
        }
    }

    double d_f = (f_max - f_min) / static_cast<double>(num_freqs - 1);
    std::vector<double> frequencies(num_freqs);
    for(size_t i = 0; i < num_freqs; i++) {
        frequencies[i] = f_min + static_cast<double>(i) * d_f;
    }
    std::vector<std::vector<double>> live_spectra(spectrometer_set.size());
    if(the_gui.running_live() && num_freqs > 0) {
        for(size_t i = 0; i < spectrometer_set.size(); i++) {
            live_spectra[i] = (*spectrometer_set[i]).acquire_single_frequency(frequencies);
        }
    } else
        return;

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }

    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = frequencies.data();
    s_plot.hasDataX = true;
    s_plot.data = live_spectra[0].data();
    s_plot.Npts = num_freqs;
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_ymin;
    s_plot.forceYmax = force_ymax;
    if(force_ymax)
        s_plot.forcedYmax = y_max;
    if(force_ymin)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Frequency (THz)";
    s_plot.yLabel = "Spectrum (counts)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    s_plot.ExtraLines = static_cast<int>(live_spectra.size()) - 1;
    if(live_spectra.size() > 1)
        s_plot.data2 = live_spectra[1].data();
    if(live_spectra.size() > 2)
        s_plot.data3 = live_spectra[2].data();
    if(live_spectra.size() > 3)
        s_plot.data4 = live_spectra[3].data();
    s_plot.plot(cr);
}

void draw_interference_spectrum(GtkDrawingArea *area,
                                cairo_t *cr,
                                int width,
                                int height,
                                gpointer data) {
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    if(!the_interference_controller.check_configuration_status())
        handle_reset_controller();
    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_y = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_y = true;
    }

    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_InterferenceF.svg");
        s_plot.SVGPath = svg_path;
    }

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }

    if(the_gui.running_live() &&
       !(*spectrometer_set[the_gui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[the_gui.pulldowns[0].get_value()])
            .set_integration_time(
                (unsigned long)round(1000 * the_gui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(the_gui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_interferogram(
            (*spectrometer_set[the_gui.pulldowns[0].get_value()]));
    }
    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = the_interference_controller.get_frequencies();
    s_plot.hasDataX = true;
    s_plot.data = the_interference_controller.get_interference_data();
    s_plot.Npts = the_interference_controller.get_num_freqs();
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_y;
    s_plot.forceYmax = force_y;
    s_plot.forcedYmax = y_max;
    if(force_y)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Frequency (THz)";
    s_plot.yLabel = "Spectrum (counts)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    s_plot.ExtraLines = 0;
    if(the_interference_controller.get_reference_a_vector().size() ==
       the_interference_controller.get_num_freqs()) {
        s_plot.data2 = the_interference_controller.get_reference_a();
        s_plot.ExtraLines = 1;
        if(the_interference_controller.get_reference_b_vector().size() ==
           the_interference_controller.get_num_freqs()) {
            s_plot.data3 = the_interference_controller.get_reference_b();
            s_plot.ExtraLines = 2;
        }
    }
    s_plot.plot(cr);
}

void draw_interference_spectrum_time(GtkDrawingArea *area,
                                     cairo_t *cr,
                                     int width,
                                     int height,
                                     gpointer data) {
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_y = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_y = true;
    }

    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_InterferenceF.svg");
        s_plot.SVGPath = svg_path;
    }

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }
    double t0 = the_gui.text_boxes[18].value_double();
    double sigma = the_gui.text_boxes[19].value_double();
    double ord = the_gui.text_boxes[20].value_double();
    the_interference_controller.set_time_filter(t0, sigma, ord);
    if(the_gui.running_live() &&
       !(*spectrometer_set[the_gui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[the_gui.pulldowns[0].get_value()])
            .set_integration_time(
                (unsigned long)round(1000 * the_gui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(the_gui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_interferogram(
            (*spectrometer_set[the_gui.pulldowns[0].get_value()]));
    }
    the_interference_controller.generate_time_plot();

    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = the_interference_controller.get_time_scale();
    s_plot.hasDataX = true;
    s_plot.data = the_interference_controller.get_time_data();
    s_plot.Npts = the_interference_controller.get_num_times();
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_y;
    s_plot.forceYmax = force_y;
    s_plot.forcedYmax = y_max;
    if(force_y)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Time (fs)";
    s_plot.yLabel = "Spectrum (counts)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    s_plot.data2 = the_interference_controller.get_time_reference_a();
    s_plot.data3 = the_interference_controller.get_time_reference_b();
    s_plot.data4 = the_interference_controller.get_time_filter();
    s_plot.ExtraLines = 3;
    s_plot.plot(cr);
}

void draw_interference_phase(GtkDrawingArea *area,
                             cairo_t *cr,
                             int width,
                             int height,
                             gpointer data) {
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_y = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_y = true;
    }

    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_InterferencePhase.svg");
        s_plot.SVGPath = svg_path;
    }

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }

    double t0 = the_gui.text_boxes[18].value_double();
    double sigma = the_gui.text_boxes[19].value_double();
    double ord = the_gui.text_boxes[20].value_double();
    the_interference_controller.set_time_filter(t0, sigma, ord);
    if(the_gui.running_live() &&
       !(*spectrometer_set[the_gui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[the_gui.pulldowns[0].get_value()])
            .set_integration_time(
                (unsigned long)round(1000 * the_gui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(the_gui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_phase(
            (*spectrometer_set[the_gui.pulldowns[0].get_value()]));
    }

    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = the_interference_controller.get_frequencies();
    s_plot.hasDataX = true;
    s_plot.data = the_interference_controller.get_mean_phase_data();
    s_plot.Npts = the_interference_controller.get_num_freqs();
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_y;
    s_plot.forceYmax = force_y;
    s_plot.forcedYmax = y_max;
    if(force_y)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Frequency (THz)";
    s_plot.yLabel = "Phase (rad)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    s_plot.ExtraLines = 0;
    s_plot.plot(cr);
}

void draw_interference_group_delay(GtkDrawingArea *area,
                                   cairo_t *cr,
                                   int width,
                                   int height,
                                   gpointer data) {
    int active_spectrometer = the_gui.pulldowns[0].get_value();
    if(active_spectrometer < 0)
        return;
    handle_reset_controller();
    LwePlot s_plot;
    bool save_svg = the_gui.save_svg > 0;
    if(save_svg) {
        the_gui.save_svg--;
    }
    bool log_plot = false;
    if(the_gui.checkboxes[1].is_checked()) {
        log_plot = true;
    }

    bool force_x = false;
    double x_min = the_gui.text_boxes[22].value_double();
    double x_max = the_gui.text_boxes[23].value_double();
    if(x_min != x_max && x_max > x_min) {
        force_x = true;
    }
    bool force_y = false;
    double y_min = the_gui.text_boxes[24].value_double();
    double y_max = the_gui.text_boxes[25].value_double();
    if(y_min != y_max && y_max > y_min) {
        force_y = true;
    }

    if(save_svg) {
        std::string svg_path;
        the_gui.file_paths[0].copy_buffer(svg_path);
        svg_path.append("_InterferencePhase.svg");
        s_plot.SVGPath = svg_path;
    }

    LweColor main_color(0.5, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[1].value_double() + the_gui.text_boxes[2].value_double() +
               the_gui.text_boxes[3].value_double())) {
        main_color = LweColor(the_gui.text_boxes[1].value_double(),
                              the_gui.text_boxes[2].value_double(),
                              the_gui.text_boxes[3].value_double(),
                              1);
    }

    LweColor over_lay0_color(1, 0, 1, 1);
    if(0.0 != (the_gui.text_boxes[4].value_double() + the_gui.text_boxes[5].value_double() +
               the_gui.text_boxes[6].value_double())) {
        over_lay0_color = LweColor(the_gui.text_boxes[4].value_double(),
                                   the_gui.text_boxes[5].value_double(),
                                   the_gui.text_boxes[6].value_double(),
                                   1);
    }

    LweColor over_lay1_color(1, 0.5, 0, 1);
    if(0.0 != (the_gui.text_boxes[7].value_double() + the_gui.text_boxes[8].value_double() +
               the_gui.text_boxes[9].value_double())) {
        over_lay1_color = LweColor(the_gui.text_boxes[7].value_double(),
                                   the_gui.text_boxes[8].value_double(),
                                   the_gui.text_boxes[9].value_double(),
                                   1);
    }

    LweColor over_lay2_color(0, 1, 1, 1);
    if(0.0 != (the_gui.text_boxes[10].value_double() + the_gui.text_boxes[11].value_double() +
               the_gui.text_boxes[12].value_double())) {
        over_lay2_color = LweColor(the_gui.text_boxes[10].value_double(),
                                   the_gui.text_boxes[11].value_double(),
                                   the_gui.text_boxes[12].value_double(),
                                   1);
    }

    double t0 = the_gui.text_boxes[18].value_double();
    double sigma = the_gui.text_boxes[19].value_double();
    double ord = the_gui.text_boxes[20].value_double();
    the_interference_controller.set_time_filter(t0, sigma, ord);
    if(the_gui.running_live() &&
       !(*spectrometer_set[the_gui.pulldowns[0].get_value()]).check_lock()) {
        (*spectrometer_set[the_gui.pulldowns[0].get_value()])
            .set_integration_time(
                (unsigned long)round(1000 * the_gui.text_boxes[0].value_double()));
        the_interference_controller.set_averaging(the_gui.checkboxes[3].is_checked());
        the_interference_controller.acquire_new_phase(
            (*spectrometer_set[the_gui.pulldowns[0].get_value()]));
    }

    s_plot.height = height;
    s_plot.width = width;
    s_plot.dataX = the_interference_controller.get_frequencies();
    s_plot.hasDataX = true;
    s_plot.data = the_interference_controller.get_group_delay_data();
    s_plot.Npts = the_interference_controller.get_num_freqs();
    s_plot.logScale = log_plot;
    s_plot.forceYmin = force_y;
    s_plot.forceYmax = force_y;
    s_plot.forcedYmax = y_max;
    if(force_y)
        s_plot.forcedYmin = y_min;
    s_plot.color = main_color;
    s_plot.color2 = over_lay0_color;
    s_plot.color3 = over_lay1_color;
    s_plot.color4 = over_lay2_color;
    s_plot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    s_plot.xLabel = "Frequency (THz)";
    s_plot.yLabel = "Group delay (ps)";
    s_plot.forceXmax = force_x;
    s_plot.forceXmin = force_x;
    s_plot.forcedXmax = x_max;
    s_plot.forcedXmin = x_min;
    s_plot.markers = false;
    s_plot.ExtraLines = 0;
    s_plot.plot(cr);
}

bool update_display() {
    the_gui.request_plot_update();
    the_gui.console.update_from_buffer();
    the_gui.apply_update();
    return true;
}

void save_path_callback() {
#ifdef __APPLE__
    theGui.path_buffer = pathFromAppleSaveDialog();
#endif
    the_gui.request_save_path_update();
}

static void activate(GtkApplication *app, gpointer user_data) {
#if defined __linux__ || defined __APPLE__
    setlocale(LC_NUMERIC, "en_US.UTF-8");
#else
    setlocale(LC_NUMERIC, "en_US");
#endif
    the_gui.activate(app);
}

int main(int argc, char **argv) {
    GtkApplication *app =
        gtk_application_new("io.github.NickKarpowicz.Scarab", (GApplicationFlags)0);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    int status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    return status;
}

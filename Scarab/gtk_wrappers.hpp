#pragma once
#include <algorithm>
#include <gcem.hpp>
#include <gtk/gtk.h>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>
#ifdef __APPLE__
#import <Cocoa/Cocoa.h>
#include <mach-o/dyld.h>
#endif
// GLOBAL VARIABLE: GTK MUTEX
static std::mutex gt_kmutex;
#include "ExternalLibraries/LightwaveExplorerPlots.h"
#include <format>

class GtkGuiElement {
  public:
    GtkWidget *label{};
    GtkWidget *element_handle{};
    int m_x{};
    int m_y{};
    int m_width{};
    int m_height{};
    bool is_attached = false;
    GtkWidget *grid{};

    void set_position(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        std::unique_lock gt_klock(gt_kmutex);
        if(grid)
            gtk_grid_remove(GTK_GRID(grid), element_handle);
        grid = grid;
        m_x = x;
        m_y = y;
        m_width = width;
        m_height = height;
        gtk_grid_attach(GTK_GRID(grid), element_handle, m_x, m_y, m_width, m_height);
    }
    void remove() { gtk_grid_remove(GTK_GRID(grid), element_handle); }
    void set_label(const int x, const int y, const char *label_text) {
        std::unique_lock gt_klock(gt_kmutex);
        label = gtk_label_new(label_text);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), 45);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(grid), label, m_x + x, m_y + y, 6, 1);
    }
    void set_label(const int x,
                  const int y,
                  const char *label_text,
                  const int characters,
                  const int grids) {
        std::unique_lock gt_klock(gt_kmutex);
        label = gtk_label_new(label_text);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), characters);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(grid), label, m_x + x, m_y + y, grids, 1);
    }
    void squeeze() {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_widget_set_valign(element_handle, GTK_ALIGN_START);
        gtk_widget_set_halign(element_handle, GTK_ALIGN_END);
    }
    void vertical_thick() {
        if(!GTK_IS_WIDGET(element_handle))
            return;
        std::unique_lock gt_klock(gt_kmutex);
        gtk_widget_set_valign(element_handle, GTK_ALIGN_FILL);
    }
    void set_tooltip(const char *tooltip_text) {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_widget_set_tooltip_text(element_handle, tooltip_text);
    }
};

class GtkTextBox : public GtkGuiElement {
  public:
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        std::unique_lock gt_klock(gt_kmutex);
        element_handle = gtk_entry_new();
        // gtk_widget_set_halign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_hexpand(element_handle, false);
        gtk_widget_set_vexpand(element_handle, false);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(element_handle), 8);
        gt_klock.unlock();
        set_position(grid, x, y, width, height);
    }

    double value_double() {
        std::unique_lock gt_klock(gt_kmutex);
        double sdata;
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> sdata;
            return sdata;
        } else {
            return 0.;
        }
    }

    int value_int() const {
        std::unique_lock gt_klock(gt_kmutex);
        int sdata;
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> sdata;
            return sdata;
        } else {
            return 0;
        }
    }

    void value_to_pointer(int *sdata) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void value_to_pointer(int64_t *sdata) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void value_to_pointer(double *sdata) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void value_to_two_pointers(double *sdata, double *sdata2) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        char c;
        *sdata2 = 0.0;
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata >> c >> *sdata2;
        }
    }

    void value_to_two_pointers(double multiplier, double *sdata, double *sdata2) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        char c;
        *sdata2 = 0.0;
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata >> c >> *sdata2;
            *sdata *= multiplier;
            *sdata2 *= multiplier;
        }
    }

    void value_to_pointer(double multiplier, double *sdata) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
            *sdata *= multiplier;
        }
    }

    void set_max_characters(const int char_limit) {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(element_handle), char_limit);
    }

    void set_to_double(const double in) {
        std::unique_lock gt_klock(gt_kmutex);
        std::string s = std::format(std::string_view("{:g}"), in);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }
    template <typename... Args> void overwrite_print(std::string_view format, Args &&...args) {
        std::unique_lock gt_klock(gt_kmutex);
        std::string s = std::vformat(format, std::make_format_args(args...));
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }

    void copy_buffer(char *destination, const int64_t max_length) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        std::string s(gtk_entry_buffer_get_text(buf));
        s.append("\0");
        s.copy(destination, max_length);
    }

    void copy_buffer(std::string &destination) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(element_handle));
        std::string s(gtk_entry_buffer_get_text(buf));
        destination = s;
    }
};

static gboolean scroll_text_view_to_end_handler(gpointer data) {
    std::unique_lock gt_klock(gt_kmutex);
    GtkAdjustment *adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(data));
    gtk_adjustment_set_value(adjustment, gtk_adjustment_get_upper(adjustment));
    return false;
}

class GtkConsole : public GtkGuiElement {
    GtkWidget *console_text{};
    bool has_new_text{};
    GtkTextBuffer *buf{};

  public:
    std::string text_buffer;
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        console_text = gtk_text_view_new();
        gtk_text_view_set_accepts_tab(GTK_TEXT_VIEW(console_text), false);
        element_handle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(element_handle), console_text);
        gtk_widget_set_vexpand(element_handle, true);
        gtk_widget_set_vexpand(console_text, true);
        gtk_widget_set_hexpand(element_handle, true);
        gtk_widget_set_hexpand(console_text, true);
        set_position(grid, x, y, width, height);
        buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_text));

        gtk_text_buffer_create_tag(buf, "function", "foreground", "#00FFFFFF", NULL);
        gtk_text_buffer_create_tag(buf, "comment", "foreground", "#006600FF", NULL);
        gtk_text_buffer_create_tag(buf, "variable", "foreground", "#FF00FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "delegate", "foreground", "#FF8800FF", NULL);
        gtk_text_buffer_create_tag(buf, "interface", "foreground", "#FF0088FF", NULL);
        gtk_text_buffer_create_tag(buf, "error", "foreground", "#FF0000FF", NULL);
        gtk_text_buffer_create_tag(buf, "parenthesis", "foreground", "#CC99FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "string", "foreground", "#FFAA00FF", NULL);
    }
    void init(GtkWidget *grid, int x, int y, int width, int height, int min_width, int min_height) {
        console_text = gtk_text_view_new();
        element_handle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(element_handle), min_height);
        gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(element_handle), min_width);
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(element_handle), console_text);
        set_position(grid, x, y, width, height);
        buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_text));
        gtk_text_buffer_create_tag(buf, "function", "foreground", "#00FFFFFF", NULL);
        gtk_text_buffer_create_tag(buf, "comment", "foreground", "#006600FF", NULL);
        gtk_text_buffer_create_tag(buf, "variable", "foreground", "#FF00FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "delegate", "foreground", "#FF8800FF", NULL);
        gtk_text_buffer_create_tag(buf, "interface", "foreground", "#FF0088FF", NULL);
        gtk_text_buffer_create_tag(buf, "error", "foreground", "#FF0000FF", NULL);
        gtk_text_buffer_create_tag(buf, "parenthesis", "foreground", "#9900CCFF", NULL);
    }

    template <typename... Args> void c_print(std::string_view format, Args &&...args) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        std::string s = std::vformat(format, std::make_format_args(args...));
        gtk_text_buffer_insert_markup(buf, &stop, s.c_str(), -1);
        gt_klock.unlock();
        scroll_to_end();
    }

    template <typename... Args> void t_print(std::string_view format, Args &&...args) {
        std::unique_lock gt_klock(gt_kmutex);
        std::string s = std::vformat(format, std::make_format_args(args...));
        text_buffer.append(s);
        has_new_text = true;
    }
    void scroll_to_end() {
        std::unique_lock gt_klock(gt_kmutex);
        g_idle_add_full(G_PRIORITY_DEFAULT_IDLE, scroll_text_view_to_end_handler, element_handle, NULL);
    }
    void update_from_buffer() {
        if(has_new_text) {
            std::unique_lock gt_klock(gt_kmutex);
            has_new_text = false;
            GtkTextIter end;
            gtk_text_buffer_get_end_iter(buf, &end);
            gtk_text_buffer_insert_markup(buf, &end, text_buffer.c_str(), -1);
            text_buffer.clear();
            gt_klock.unlock();
            scroll_to_end();
        }
    }

    void direct_overwrite_print(const char *s_in) {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_text_buffer_set_text(buf, s_in, -1);
        text_buffer.clear();
        gt_klock.unlock();
        scroll_to_end();
    }

    template <typename... Args> void overwrite_print(std::string_view format, Args &&...args) {
        std::unique_lock gt_klock(gt_kmutex);
        std::string s = Svformat(format, Smake_format_args(args...));
        gtk_text_buffer_set_text(buf, s.c_str(), (int)s.length());
        text_buffer.clear();
        gt_klock.unlock();
        scroll_to_end();
    }

    void copy_buffer(char *destination, int64_t max_length) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        char *real_buf = gtk_text_buffer_get_text(buf, &start, &stop, false);
        std::string s(real_buf);
        s.copy(destination, max_length);
    }

    void copy_buffer(std::string &s) {
        std::unique_lock gt_klock(gt_kmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        char *real_buf = gtk_text_buffer_get_text(buf, &start, &stop, false);
        std::string c(real_buf);
        s = c;
    }

    void clear() {
        std::unique_lock gt_klock(gt_kmutex);
        char empty_buffer[] = "";
        text_buffer.clear();
        GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_text));
        gtk_text_buffer_set_text(buf, empty_buffer, 0);
        gt_klock.unlock();
        scroll_to_end();
    }
};

class GtkPushbutton : public GtkGuiElement {
  public:
    void init(const char *button_name,
              GtkWidget *grid,
              int x,
              int y,
              int width,
              int height,
              auto button_function) {
        std::unique_lock gt_klock(gt_kmutex);
        element_handle = gtk_button_new_with_label(button_name);
        gtk_widget_set_hexpand(element_handle, false);
        gtk_widget_set_vexpand(element_handle, false);
        gtk_widget_set_valign(element_handle, GTK_ALIGN_START);
        gt_klock.unlock();
        set_position(grid, x, y, width, height);
        set_function(button_function);
    }
    void init(const char *button_name,
              GtkWidget *grid,
              int x,
              int y,
              int width,
              int height,
              auto button_function,
              gpointer function_data) {
        std::unique_lock gt_klock(gt_kmutex);
        element_handle = gtk_button_new_with_label(button_name);
        gtk_widget_set_hexpand(element_handle, false);
        gtk_widget_set_vexpand(element_handle, false);
        gt_klock.unlock();
        set_position(grid, x, y, width, height);
        set_function(button_function, function_data);
    }
    void set_function(auto button_function) {
        std::unique_lock gt_klock(gt_kmutex);
        g_signal_connect(element_handle, "clicked", G_CALLBACK(button_function), NULL);
    }
    void set_function(auto button_function, gpointer param) {
        std::unique_lock gt_klock(gt_kmutex);
        g_signal_connect(element_handle, "clicked", G_CALLBACK(button_function), param);
    }
};

class GtkCheck : public GtkGuiElement {
  public:
    void init(const char *button_name, GtkWidget *grid, int x, int y, int width, int height) {
        std::unique_lock gt_klock(gt_kmutex);
        element_handle = gtk_check_button_new_with_label(button_name);
        gt_klock.unlock();
        set_position(grid, x, y, width, height);
    }
    bool is_checked() {
        std::unique_lock gt_klock(gt_kmutex);
        return (bool)gtk_check_button_get_active(GTK_CHECK_BUTTON(element_handle));
    }
    void set_function(auto button_function) {
        std::unique_lock gt_klock(gt_kmutex);
        g_signal_connect_after(element_handle, "toggled", G_CALLBACK(button_function), NULL);
    }
};

class GtkPulldown : public GtkGuiElement {
    std::vector<std::string> entry_names;

  public:
    void add_element(const char *newelement) {
        std::string s(newelement);
        strip_line_breaks(s);
        entry_names.push_back(s);
    }
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        // make an array of pointers to c-strings for GTK
        std::vector<const char *> string_pointers_for_gtk;
        string_pointers_for_gtk.reserve(entry_names.size() + 1);
        std::transform(entry_names.begin(),
                       entry_names.end(),
                       std::back_inserter(string_pointers_for_gtk),
                       [](const std::string &s) { return s.c_str(); });
        string_pointers_for_gtk.push_back(nullptr);

        std::unique_lock gt_klock(gt_kmutex);
        element_handle = gtk_drop_down_new_from_strings(string_pointers_for_gtk.data());
        gtk_widget_set_hexpand(element_handle, false);
        gtk_widget_set_vexpand(element_handle, false);
        gtk_widget_set_valign(element_handle, GTK_ALIGN_START);
        gt_klock.unlock();
        set_position(grid, x, y, width, height);
    }
    int get_value() {
        std::unique_lock gt_klock(gt_kmutex);
        return (int)gtk_drop_down_get_selected(GTK_DROP_DOWN(element_handle));
    }
    void set_value(int target) {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_drop_down_set_selected(GTK_DROP_DOWN(element_handle), target);
    }
    inline void remove_character_from_string(std::string &s, const char removed_char) {
        std::erase(s, removed_char);
    }
    void strip_line_breaks(std::string &s) {
        remove_character_from_string(s, '\r');
        remove_character_from_string(s, '\n');
    }
};

class GtkMainWindow {
    GtkWidget *grid = nullptr;
    GtkWidget *big_grid = nullptr;
    GtkWidget *console_grid = nullptr;
    GtkWidget *plot_grid = nullptr;
    GtkWidget *plot_controls_grid = nullptr;
    GtkWidget *console_controls_grid = nullptr;
    GtkWidget *console_controls_subgrid1 = nullptr;
    GtkWidget *console_controls_subgrid2 = nullptr;
    GtkWidget *plot_controls_subgrid1 = nullptr;
    GtkWidget *plot_controls_subgrid2 = nullptr;
    unsigned int updater_id = 0;
    bool showing_controls_panel = true;

  public:
    GtkWidget *window = nullptr;
    void init(GtkApplication *app_handle, const char *window_name, int width, int height) {
        std::unique_lock gt_klock(gt_kmutex);

        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", true, NULL);
        GtkCssProvider *text_provider = gtk_css_provider_new();
// override style more aggressively if it will be Adwaita
#if defined _WIN32 || defined __APPLE__ || defined LWEFLATPAK
        std::string styleString("label, scale { font-family: Arial; font-weight: bold; }\n "
                                "button, entry, textview { font-family: Arial; font-weight: bold; "
                                "color: #FFFFFF; background-color: #151515; }");
#else
        std::string style_string(
            "label, scale { font-family: Arial; font-weight: bold; }\n "
            "button, entry, textview { font-family: Arial; font-weight: bold;}");
#endif

#if defined __APPLE__
        gtk_css_provider_load_from_data(textProvider, styleString.c_str(), -1);
#else
        gtk_css_provider_load_from_string(text_provider, style_string.c_str());
#endif
        gtk_style_context_add_provider_for_display(gdk_display_get_default(),
                                                   GTK_STYLE_PROVIDER(text_provider),
                                                   GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        GtkCssProvider *button_shrinker = gtk_css_provider_new();
        std::string button_style("label, scale, button, entry, textview "
                                "{ min-height: 17px; min-width: 8px; }");
#if defined __APPLE__
        gtk_css_provider_load_from_data(buttonShrinker, buttonStyle.c_str(), -1);
#else
        gtk_css_provider_load_from_string(button_shrinker, button_style.c_str());
#endif
        gtk_style_context_add_provider_for_display(gdk_display_get_default(),
                                                   GTK_STYLE_PROVIDER(button_shrinker),
                                                   GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        window = gtk_application_window_new(app_handle);
        gtk_window_set_title(GTK_WINDOW(window), window_name);
        gtk_window_set_default_size(GTK_WINDOW(window), width, height);
        big_grid = gtk_grid_new();
        console_grid = gtk_grid_new();
        console_controls_grid = gtk_grid_new();
        console_controls_subgrid1 = gtk_grid_new();
        console_controls_subgrid2 = gtk_grid_new();
        plot_grid = gtk_grid_new();
        plot_controls_grid = gtk_grid_new();
        plot_controls_subgrid1 = gtk_grid_new();
        plot_controls_subgrid2 = gtk_grid_new();
        grid = gtk_grid_new();
        gtk_grid_set_row_spacing(GTK_GRID(grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(grid), 1);

        gtk_grid_set_row_spacing(GTK_GRID(big_grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(big_grid), 1);
        gtk_grid_set_row_spacing(GTK_GRID(plot_grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(plot_grid), 1);
        gtk_grid_set_row_homogeneous(GTK_GRID(grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(console_controls_grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(plot_grid), true);
        gtk_grid_set_row_homogeneous(GTK_GRID(plot_grid), true);
        gtk_grid_set_row_spacing(GTK_GRID(console_grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(console_controls_grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(console_grid), 1);

        gtk_window_set_child(GTK_WINDOW(window), big_grid);

        gtk_widget_set_hexpand(grid, false);
        gtk_widget_set_vexpand(grid, false);
        gtk_widget_set_halign(grid, GTK_ALIGN_END);
        gtk_widget_set_valign(grid, GTK_ALIGN_START);

        gtk_widget_set_hexpand(console_grid, false);
        gtk_widget_set_vexpand(console_grid, true);
        gtk_widget_set_halign(console_grid, GTK_ALIGN_FILL);
        gtk_widget_set_valign(console_grid, GTK_ALIGN_FILL);

        gtk_widget_set_hexpand(console_controls_grid, false);
        gtk_widget_set_halign(console_controls_grid, GTK_ALIGN_FILL);

        gtk_widget_set_valign(console_controls_subgrid1, GTK_ALIGN_CENTER);
        gtk_widget_set_halign(console_controls_subgrid2, GTK_ALIGN_END);
        gtk_grid_attach(GTK_GRID(big_grid), console_grid, 0, 1, 1, 1);
        gtk_grid_attach(GTK_GRID(big_grid), console_controls_grid, 0, 2, 1, 1);
        gtk_grid_attach(GTK_GRID(console_controls_grid), console_controls_subgrid1, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(console_controls_grid), console_controls_subgrid2, 1, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(big_grid), plot_grid, 1, 0, 1, 2);
        gtk_grid_attach(GTK_GRID(big_grid), grid, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(big_grid), plot_controls_grid, 1, 2, 1, 1);
        gtk_grid_attach(GTK_GRID(plot_controls_grid), plot_controls_subgrid1, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(plot_controls_grid), plot_controls_subgrid2, 1, 0, 1, 1);
    }
    void present() {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_window_present(GTK_WINDOW(window));
    }

    GtkWidget *parent_handle() const { return grid; }
    void toggle_settings_panel() {
        std::unique_lock gt_klock(gt_kmutex);
        showing_controls_panel = !showing_controls_panel;
        gtk_widget_set_visible(grid, showing_controls_panel);
        gtk_widget_set_visible(console_grid, showing_controls_panel);
        gtk_widget_set_visible(console_controls_subgrid1, showing_controls_panel);
        gtk_widget_set_visible(console_controls_subgrid2, showing_controls_panel);
    }
    GtkWidget *parent_handle(const int index) const {
        switch(index) {
            case 0:
                return grid;
            case 1:
                return console_grid;
            case 2:
                return plot_grid;
            case 3:
                return plot_controls_subgrid1;
            case 4:
                return plot_controls_subgrid2;
            case 5:
                return console_controls_subgrid1;
            case 6:
                return console_controls_subgrid2;
        }
        return grid;
    }

    GtkWindow *window_handle() { return GTK_WINDOW(window); }
    void connect_update_function(auto function) {
        updater_id = g_timeout_add(50, G_SOURCE_FUNC(function), NULL);
    }
    void remove_update_function() { g_source_remove(updater_id); }
};

class GtkDrawBox : public GtkGuiElement {
  public:
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        std::unique_lock gt_klock(gt_kmutex);
        element_handle = gtk_drawing_area_new();
        gtk_widget_set_hexpand(element_handle, true);
        gtk_widget_set_vexpand(element_handle, true);
        gt_klock.unlock();
        set_position(grid, x, y, width, height);
        gt_klock.lock();
        gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(element_handle), 320);
        gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(element_handle), 12);
    }
    void set_drawing_function(GtkDrawingAreaDrawFunc the_function) {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(element_handle), the_function, NULL, NULL);
    }
    void queue_draw() {
        if(!GTK_IS_WIDGET(element_handle))
            return;
        std::unique_lock gt_klock(gt_kmutex);
        gtk_widget_queue_draw(element_handle);
    }
    void no_vertical_expantion() {
        std::unique_lock gt_klock(gt_kmutex);
        gtk_widget_set_vexpand(element_handle, false);
    }
};

static void path_from_load_dialog_callback(GObject *gobject, GAsyncResult *result, gpointer data) {
    GtkTextBox &destination_path_box = *reinterpret_cast<GtkTextBox *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        std::string result_string(g_file_get_path(file));
        destination_path_box.overwrite_print(result_string);
    }
}
static void
path_from_load_dialog_to_string_callback(GObject *gobject, GAsyncResult *result, gpointer data) {
    std::string &destination_path = *reinterpret_cast<std::string *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(gobject), result, &error);

    if(error == nullptr) {
        destination_path = std::string(g_file_get_path(file));
    } else {
        destination_path = std::string("?LWE_NOPATH??");
    }
}

inline void path_from_load_dialog(GtkTextBox &destination_path_box) {
    GtkFileDialog *dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(dialog, NULL, NULL, path_from_load_dialog_callback, &destination_path_box);
}
inline void path_from_load_dialog(std::string &destination_path) {
    GtkFileDialog *dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(dialog, NULL, NULL, path_from_load_dialog_to_string_callback, &destination_path);
}

[[maybe_unused]] static void
path_from_save_dialog_callback(GObject *gobject, GAsyncResult *result, gpointer data) {
    GtkTextBox &destination_path_box = *reinterpret_cast<GtkTextBox *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        std::string path(g_file_get_path(file));
        destination_path_box.overwrite_print(path);
    }
}
[[maybe_unused]] static void
path_from_save_dialog_string_callback(GObject *gobject, GAsyncResult *result, gpointer data) {
    std::string &destination_path = *reinterpret_cast<std::string *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        destination_path = std::string(g_file_get_path(file));
    }
}
inline void path_from_save_dialog(GtkTextBox &destination_path_box) {
#ifdef __APPLE__
    NSString *filePath;
    NSSavePanel *savePanel = [NSSavePanel savePanel];
    if([savePanel runModal] == NSModalResponseOK) {
        filePath = [savePanel URL].path;
        destinationPathBox.overwritePrint("{}", [filePath UTF8String]);
    }
#else
    GtkFileDialog *dialog = gtk_file_dialog_new();
    GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
    GtkFileFilter *filter = gtk_file_filter_new();

    gtk_file_filter_add_suffix(filter, "zip");
    gtk_file_filter_set_name(filter, "Compressed (.zip)");
    g_list_store_append(filters, filter);

    filter = gtk_file_filter_new();
    gtk_file_filter_add_pattern(filter, "*");
    gtk_file_filter_set_name(filter, "All Files");
    g_list_store_append(filters, filter);

    gtk_file_dialog_set_filters(dialog, G_LIST_MODEL(filters));
    gtk_file_dialog_save(dialog, NULL, NULL, path_from_save_dialog_callback, &destination_path_box);
#endif
}
#ifdef __APPLE__
std::string pathFromAppleSaveDialog() {
    NSString *filePath;
    NSSavePanel *savePanel = [NSSavePanel savePanel];
    if([savePanel runModal] == NSModalResponseOK) {
        filePath = [savePanel URL].path;
        return std::string([filePath UTF8String]);
    } else {
        return std::string("?LWE_NOPATH??");
    }
}
#endif
inline void path_from_save_dialog(std::string &destination_path,
                               const std::string &suffix,
                               const std::string &filetype_name) {
#ifdef __APPLE__
    NSString *filePath;
    NSSavePanel *savePanel = [NSSavePanel savePanel];
    if([savePanel runModal] == NSModalResponseOK) {
        filePath = [savePanel URL].path;
        destinationPath = [filePath UTF8String];
    }
#else
    GtkFileDialog *dialog = gtk_file_dialog_new();
    GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);

    GtkFileFilter *filter = gtk_file_filter_new();
    gtk_file_filter_add_suffix(filter, suffix.c_str());
    gtk_file_filter_set_name(filter, filetype_name.c_str());
    g_list_store_append(filters, filter);

    filter = gtk_file_filter_new();
    gtk_file_filter_add_pattern(filter, "*");
    gtk_file_filter_set_name(filter, "All Files");
    g_list_store_append(filters, filter);

    gtk_file_dialog_set_filters(dialog, G_LIST_MODEL(filters));
    gtk_file_dialog_save(dialog, NULL, NULL, path_from_save_dialog_string_callback, &destination_path);
#endif
}

typedef void (*loadingFunction)(std::string);
static void load_data_callback(GObject *gobject, GAsyncResult *result, gpointer data) {
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        std::string path(g_file_get_path(file));

        loadingFunction loading_function_pointer = reinterpret_cast<loadingFunction>(data);
        (loading_function_pointer)(path);
    }
}

inline void load_from_load_dialog(loadingFunction loading_function_pointer) {
    GtkFileDialog *dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(dialog,
                         NULL,
                         NULL,
                         load_data_callback,
                         reinterpret_cast<gpointer>(loading_function_pointer));
}

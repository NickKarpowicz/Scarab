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
static std::mutex GTKmutex;
#include "ExternalLibraries/LightwaveExplorerPlots.h"
#include <format>

class GtkGuiElement {
  public:
    GtkWidget *label{};
    GtkWidget *elementHandle{};
    int m_x{};
    int m_y{};
    int m_width{};
    int m_height{};
    bool isAttached = false;
    GtkWidget *_grid{};

    void setPosition(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        std::unique_lock GTKlock(GTKmutex);
        if(_grid)
            gtk_grid_remove(GTK_GRID(_grid), elementHandle);
        _grid = grid;
        m_x = x;
        m_y = y;
        m_width = width;
        m_height = height;
        gtk_grid_attach(GTK_GRID(_grid), elementHandle, m_x, m_y, m_width, m_height);
    }
    void remove() { gtk_grid_remove(GTK_GRID(_grid), elementHandle); }
    void setLabel(const int x, const int y, const char *labelText) {
        std::unique_lock GTKlock(GTKmutex);
        label = gtk_label_new(labelText);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), 45);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(_grid), label, m_x + x, m_y + y, 6, 1);
    }
    void setLabel(const int x,
                  const int y,
                  const char *labelText,
                  const int characters,
                  const int grids) {
        std::unique_lock GTKlock(GTKmutex);
        label = gtk_label_new(labelText);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), characters);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(_grid), label, m_x + x, m_y + y, grids, 1);
    }
    void squeeze() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_halign(elementHandle, GTK_ALIGN_END);
    }
    void vertical_thick() {
        if(!GTK_IS_WIDGET(elementHandle))
            return;
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_FILL);
    }
    void setTooltip(const char *tooltipText) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_tooltip_text(elementHandle, tooltipText);
    }
};

class GtkTextBox : public GtkGuiElement {
  public:
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_entry_new();
        // gtk_widget_set_halign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(elementHandle), 8);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }

    double valueDouble() {
        std::unique_lock GTKlock(GTKmutex);
        double sdata;
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
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

    int valueInt() const {
        std::unique_lock GTKlock(GTKmutex);
        int sdata;
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
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

    void valueToPointer(int *sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToPointer(int64_t *sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToPointer(double *sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToTwoPointers(double *sdata, double *sdata2) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        char c;
        *sdata2 = 0.0;
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata >> c >> *sdata2;
        }
    }

    void valueToTwoPointers(double multiplier, double *sdata, double *sdata2) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
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

    void valueToPointer(double multiplier, double *sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if(len > 0) {
            char *cbuf = (char *)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
            *sdata *= multiplier;
        }
    }

    void setMaxCharacters(const int charLimit) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(elementHandle), charLimit);
    }

    void setToDouble(const double in) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = std::format(std::string_view("{:g}"), in);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }
    template <typename... Args> void overwritePrint(std::string_view format, Args &&...args) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = std::vformat(format, std::make_format_args(args...));
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }

    void copyBuffer(char *destination, const int64_t maxLength) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        std::string s(gtk_entry_buffer_get_text(buf));
        s.append("\0");
        s.copy(destination, maxLength);
    }

    void copyBuffer(std::string &destination) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer *buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        std::string s(gtk_entry_buffer_get_text(buf));
        destination = s;
    }
};

static gboolean scrollTextViewToEndHandler(gpointer data) {
    std::unique_lock GTKlock(GTKmutex);
    GtkAdjustment *adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(data));
    gtk_adjustment_set_value(adjustment, gtk_adjustment_get_upper(adjustment));
    return false;
}

class GtkConsole : public GtkGuiElement {
    GtkWidget *consoleText{};
    bool hasNewText{};
    GtkTextBuffer *buf{};

  public:
    std::string textBuffer;
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        consoleText = gtk_text_view_new();
        gtk_text_view_set_accepts_tab(GTK_TEXT_VIEW(consoleText), false);
        elementHandle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(elementHandle), consoleText);
        gtk_widget_set_vexpand(elementHandle, true);
        gtk_widget_set_vexpand(consoleText, true);
        gtk_widget_set_hexpand(elementHandle, true);
        gtk_widget_set_hexpand(consoleText, true);
        setPosition(grid, x, y, width, height);
        buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));

        gtk_text_buffer_create_tag(buf, "function", "foreground", "#00FFFFFF", NULL);
        gtk_text_buffer_create_tag(buf, "comment", "foreground", "#006600FF", NULL);
        gtk_text_buffer_create_tag(buf, "variable", "foreground", "#FF00FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "delegate", "foreground", "#FF8800FF", NULL);
        gtk_text_buffer_create_tag(buf, "interface", "foreground", "#FF0088FF", NULL);
        gtk_text_buffer_create_tag(buf, "error", "foreground", "#FF0000FF", NULL);
        gtk_text_buffer_create_tag(buf, "parenthesis", "foreground", "#CC99FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "string", "foreground", "#FFAA00FF", NULL);
    }
    void init(GtkWidget *grid, int x, int y, int width, int height, int minWidth, int minHeight) {
        consoleText = gtk_text_view_new();
        elementHandle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(elementHandle), minHeight);
        gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(elementHandle), minWidth);
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(elementHandle), consoleText);
        setPosition(grid, x, y, width, height);
        buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_create_tag(buf, "function", "foreground", "#00FFFFFF", NULL);
        gtk_text_buffer_create_tag(buf, "comment", "foreground", "#006600FF", NULL);
        gtk_text_buffer_create_tag(buf, "variable", "foreground", "#FF00FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "delegate", "foreground", "#FF8800FF", NULL);
        gtk_text_buffer_create_tag(buf, "interface", "foreground", "#FF0088FF", NULL);
        gtk_text_buffer_create_tag(buf, "error", "foreground", "#FF0000FF", NULL);
        gtk_text_buffer_create_tag(buf, "parenthesis", "foreground", "#9900CCFF", NULL);
    }

    template <typename... Args> void cPrint(std::string_view format, Args &&...args) {
        std::unique_lock GTKlock(GTKmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        std::string s = std::vformat(format, std::make_format_args(args...));
        gtk_text_buffer_insert_markup(buf, &stop, s.c_str(), -1);
        GTKlock.unlock();
        scrollToEnd();
    }

    template <typename... Args> void tPrint(std::string_view format, Args &&...args) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = std::vformat(format, std::make_format_args(args...));
        textBuffer.append(s);
        hasNewText = true;
    }
    void scrollToEnd() {
        std::unique_lock GTKlock(GTKmutex);
        g_idle_add_full(G_PRIORITY_DEFAULT_IDLE, scrollTextViewToEndHandler, elementHandle, NULL);
    }
    void updateFromBuffer() {
        if(hasNewText) {
            std::unique_lock GTKlock(GTKmutex);
            hasNewText = false;
            GtkTextIter end;
            gtk_text_buffer_get_end_iter(buf, &end);
            gtk_text_buffer_insert_markup(buf, &end, textBuffer.c_str(), -1);
            textBuffer.clear();
            GTKlock.unlock();
            scrollToEnd();
        }
    }

    void directOverwritePrint(const char *sIn) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_text_buffer_set_text(buf, sIn, -1);
        textBuffer.clear();
        GTKlock.unlock();
        scrollToEnd();
    }

    template <typename... Args> void overwritePrint(std::string_view format, Args &&...args) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = Svformat(format, Smake_format_args(args...));
        gtk_text_buffer_set_text(buf, s.c_str(), (int)s.length());
        textBuffer.clear();
        GTKlock.unlock();
        scrollToEnd();
    }

    void copyBuffer(char *destination, int64_t maxLength) {
        std::unique_lock GTKlock(GTKmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        char *realBuf = gtk_text_buffer_get_text(buf, &start, &stop, false);
        std::string s(realBuf);
        s.copy(destination, maxLength);
    }

    void copyBuffer(std::string &s) {
        std::unique_lock GTKlock(GTKmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        char *realBuf = gtk_text_buffer_get_text(buf, &start, &stop, false);
        std::string c(realBuf);
        s = c;
    }

    void clear() {
        std::unique_lock GTKlock(GTKmutex);
        char emptyBuffer[] = "";
        textBuffer.clear();
        GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_set_text(buf, emptyBuffer, 0);
        GTKlock.unlock();
        scrollToEnd();
    }
};

class GtkPushbutton : public GtkGuiElement {
  public:
    void init(const char *buttonName,
              GtkWidget *grid,
              int x,
              int y,
              int width,
              int height,
              auto buttonFunction) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_button_new_with_label(buttonName);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
        setFunction(buttonFunction);
    }
    void init(const char *buttonName,
              GtkWidget *grid,
              int x,
              int y,
              int width,
              int height,
              auto buttonFunction,
              gpointer functionData) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_button_new_with_label(buttonName);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
        setFunction(buttonFunction, functionData);
    }
    void setFunction(auto buttonFunction) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect(elementHandle, "clicked", G_CALLBACK(buttonFunction), NULL);
    }
    void setFunction(auto buttonFunction, gpointer param) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect(elementHandle, "clicked", G_CALLBACK(buttonFunction), param);
    }
};

class GtkCheck : public GtkGuiElement {
  public:
    void init(const char *buttonName, GtkWidget *grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_check_button_new_with_label(buttonName);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
    bool isChecked() {
        std::unique_lock GTKlock(GTKmutex);
        return (bool)gtk_check_button_get_active(GTK_CHECK_BUTTON(elementHandle));
    }
    void setFunction(auto buttonFunction) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect_after(elementHandle, "toggled", G_CALLBACK(buttonFunction), NULL);
    }
};

class GtkPulldown : public GtkGuiElement {
    std::vector<std::string> entryNames;

  public:
    void addElement(const char *newelement) {
        std::string s(newelement);
        stripLineBreaks(s);
        entryNames.push_back(s);
    }
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        // make an array of pointers to c-strings for GTK
        std::vector<const char *> stringPointersForGTK;
        stringPointersForGTK.reserve(entryNames.size() + 1);
        std::transform(entryNames.begin(),
                       entryNames.end(),
                       std::back_inserter(stringPointersForGTK),
                       [](const std::string &s) { return s.c_str(); });
        stringPointersForGTK.push_back(nullptr);

        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_drop_down_new_from_strings(stringPointersForGTK.data());
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
    int getValue() {
        std::unique_lock GTKlock(GTKmutex);
        return (int)gtk_drop_down_get_selected(GTK_DROP_DOWN(elementHandle));
    }
    void setValue(int target) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_drop_down_set_selected(GTK_DROP_DOWN(elementHandle), target);
    }
    inline void removeCharacterFromString(std::string &s, const char removedChar) {
        std::erase(s, removedChar);
    }
    void stripLineBreaks(std::string &s) {
        removeCharacterFromString(s, '\r');
        removeCharacterFromString(s, '\n');
    }
};

class GtkMainWindow {
    GtkWidget *grid = nullptr;
    GtkWidget *bigGrid = nullptr;
    GtkWidget *consoleGrid = nullptr;
    GtkWidget *plotGrid = nullptr;
    GtkWidget *plotControlsGrid = nullptr;
    GtkWidget *consoleControlsGrid = nullptr;
    GtkWidget *consoleControlsSubgrid1 = nullptr;
    GtkWidget *consoleControlsSubgrid2 = nullptr;
    GtkWidget *plotControlsSubgrid1 = nullptr;
    GtkWidget *plotControlsSubgrid2 = nullptr;
    unsigned int updaterID = 0;
    bool showingControlsPanel = true;

  public:
    GtkWidget *window = nullptr;
    void init(GtkApplication *appHandle, const char *windowName, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);

        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", true, NULL);
        GtkCssProvider *textProvider = gtk_css_provider_new();
// override style more aggressively if it will be Adwaita
#if defined _WIN32 || defined __APPLE__ || defined LWEFLATPAK
        std::string styleString("label, scale { font-family: Arial; font-weight: bold; }\n "
                                "button, entry, textview { font-family: Arial; font-weight: bold; "
                                "color: #FFFFFF; background-color: #151515; }");
#else
        std::string styleString(
            "label, scale { font-family: Arial; font-weight: bold; }\n "
            "button, entry, textview { font-family: Arial; font-weight: bold;}");
#endif

#if defined __APPLE__
        gtk_css_provider_load_from_data(textProvider, styleString.c_str(), -1);
#else
        gtk_css_provider_load_from_string(textProvider, styleString.c_str());
#endif
        gtk_style_context_add_provider_for_display(gdk_display_get_default(),
                                                   GTK_STYLE_PROVIDER(textProvider),
                                                   GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        GtkCssProvider *buttonShrinker = gtk_css_provider_new();
        std::string buttonStyle("label, scale, button, entry, textview "
                                "{ min-height: 17px; min-width: 8px; }");
#if defined __APPLE__
        gtk_css_provider_load_from_data(buttonShrinker, buttonStyle.c_str(), -1);
#else
        gtk_css_provider_load_from_string(buttonShrinker, buttonStyle.c_str());
#endif
        gtk_style_context_add_provider_for_display(gdk_display_get_default(),
                                                   GTK_STYLE_PROVIDER(buttonShrinker),
                                                   GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        window = gtk_application_window_new(appHandle);
        gtk_window_set_title(GTK_WINDOW(window), windowName);
        gtk_window_set_default_size(GTK_WINDOW(window), width, height);
        bigGrid = gtk_grid_new();
        consoleGrid = gtk_grid_new();
        consoleControlsGrid = gtk_grid_new();
        consoleControlsSubgrid1 = gtk_grid_new();
        consoleControlsSubgrid2 = gtk_grid_new();
        plotGrid = gtk_grid_new();
        plotControlsGrid = gtk_grid_new();
        plotControlsSubgrid1 = gtk_grid_new();
        plotControlsSubgrid2 = gtk_grid_new();
        grid = gtk_grid_new();
        gtk_grid_set_row_spacing(GTK_GRID(grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(grid), 1);

        gtk_grid_set_row_spacing(GTK_GRID(bigGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(bigGrid), 1);
        gtk_grid_set_row_spacing(GTK_GRID(plotGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(plotGrid), 1);
        gtk_grid_set_row_homogeneous(GTK_GRID(grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(consoleControlsGrid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(plotGrid), true);
        gtk_grid_set_row_homogeneous(GTK_GRID(plotGrid), true);
        gtk_grid_set_row_spacing(GTK_GRID(consoleGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(consoleControlsGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(consoleGrid), 1);

        gtk_window_set_child(GTK_WINDOW(window), bigGrid);

        gtk_widget_set_hexpand(grid, false);
        gtk_widget_set_vexpand(grid, false);
        gtk_widget_set_halign(grid, GTK_ALIGN_END);
        gtk_widget_set_valign(grid, GTK_ALIGN_START);

        gtk_widget_set_hexpand(consoleGrid, false);
        gtk_widget_set_vexpand(consoleGrid, true);
        gtk_widget_set_halign(consoleGrid, GTK_ALIGN_FILL);
        gtk_widget_set_valign(consoleGrid, GTK_ALIGN_FILL);

        gtk_widget_set_hexpand(consoleControlsGrid, false);
        gtk_widget_set_halign(consoleControlsGrid, GTK_ALIGN_FILL);

        gtk_widget_set_valign(consoleControlsSubgrid1, GTK_ALIGN_CENTER);
        gtk_widget_set_halign(consoleControlsSubgrid2, GTK_ALIGN_END);
        gtk_grid_attach(GTK_GRID(bigGrid), consoleGrid, 0, 1, 1, 1);
        gtk_grid_attach(GTK_GRID(bigGrid), consoleControlsGrid, 0, 2, 1, 1);
        gtk_grid_attach(GTK_GRID(consoleControlsGrid), consoleControlsSubgrid1, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(consoleControlsGrid), consoleControlsSubgrid2, 1, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(bigGrid), plotGrid, 1, 0, 1, 2);
        gtk_grid_attach(GTK_GRID(bigGrid), grid, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(bigGrid), plotControlsGrid, 1, 2, 1, 1);
        gtk_grid_attach(GTK_GRID(plotControlsGrid), plotControlsSubgrid1, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(plotControlsGrid), plotControlsSubgrid2, 1, 0, 1, 1);
    }
    void present() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_window_present(GTK_WINDOW(window));
    }

    GtkWidget *parentHandle() const { return grid; }
    void toggleSettingsPanel() {
        std::unique_lock GTKlock(GTKmutex);
        showingControlsPanel = !showingControlsPanel;
        gtk_widget_set_visible(grid, showingControlsPanel);
        gtk_widget_set_visible(consoleGrid, showingControlsPanel);
        gtk_widget_set_visible(consoleControlsSubgrid1, showingControlsPanel);
        gtk_widget_set_visible(consoleControlsSubgrid2, showingControlsPanel);
    }
    GtkWidget *parentHandle(const int index) const {
        switch(index) {
            case 0:
                return grid;
            case 1:
                return consoleGrid;
            case 2:
                return plotGrid;
            case 3:
                return plotControlsSubgrid1;
            case 4:
                return plotControlsSubgrid2;
            case 5:
                return consoleControlsSubgrid1;
            case 6:
                return consoleControlsSubgrid2;
        }
        return grid;
    }

    GtkWindow *windowHandle() { return GTK_WINDOW(window); }
    void connectUpdateFunction(auto function) {
        updaterID = g_timeout_add(50, G_SOURCE_FUNC(function), NULL);
    }
    void removeUpdateFunction() { g_source_remove(updaterID); }
};

class GtkDrawBox : public GtkGuiElement {
  public:
    void init(GtkWidget *grid, const int x, const int y, const int width, const int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_drawing_area_new();
        gtk_widget_set_hexpand(elementHandle, true);
        gtk_widget_set_vexpand(elementHandle, true);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
        GTKlock.lock();
        gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(elementHandle), 320);
        gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(elementHandle), 12);
    }
    void setDrawingFunction(GtkDrawingAreaDrawFunc theFunction) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(elementHandle), theFunction, NULL, NULL);
    }
    void queueDraw() {
        if(!GTK_IS_WIDGET(elementHandle))
            return;
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_queue_draw(elementHandle);
    }
    void noVerticalExpantion() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_vexpand(elementHandle, false);
    }
};

static void pathFromLoadDialogCallback(GObject *gobject, GAsyncResult *result, gpointer data) {
    GtkTextBox &destinationPathBox = *reinterpret_cast<GtkTextBox *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        std::string resultString(g_file_get_path(file));
        destinationPathBox.overwritePrint(resultString);
    }
}
static void
pathFromLoadDialogToStringCallback(GObject *gobject, GAsyncResult *result, gpointer data) {
    std::string &destinationPath = *reinterpret_cast<std::string *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(gobject), result, &error);

    if(error == nullptr) {
        destinationPath = std::string(g_file_get_path(file));
    } else {
        destinationPath = std::string("?LWE_NOPATH??");
    }
}

inline void pathFromLoadDialog(GtkTextBox &destinationPathBox) {
    GtkFileDialog *dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(dialog, NULL, NULL, pathFromLoadDialogCallback, &destinationPathBox);
}
inline void pathFromLoadDialog(std::string &destinationPath) {
    GtkFileDialog *dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(dialog, NULL, NULL, pathFromLoadDialogToStringCallback, &destinationPath);
}

[[maybe_unused]] static void
pathFromSaveDialogCallback(GObject *gobject, GAsyncResult *result, gpointer data) {
    GtkTextBox &destinationPathBox = *reinterpret_cast<GtkTextBox *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        std::string path(g_file_get_path(file));
        destinationPathBox.overwritePrint(path);
    }
}
[[maybe_unused]] static void
pathFromSaveDialogStringCallback(GObject *gobject, GAsyncResult *result, gpointer data) {
    std::string &destinationPath = *reinterpret_cast<std::string *>(data);
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        destinationPath = std::string(g_file_get_path(file));
    }
}
inline void pathFromSaveDialog(GtkTextBox &destinationPathBox) {
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
    gtk_file_dialog_save(dialog, NULL, NULL, pathFromSaveDialogCallback, &destinationPathBox);
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
inline void pathFromSaveDialog(std::string &destinationPath,
                               const std::string &suffix,
                               const std::string &filetypeName) {
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
    gtk_file_filter_set_name(filter, filetypeName.c_str());
    g_list_store_append(filters, filter);

    filter = gtk_file_filter_new();
    gtk_file_filter_add_pattern(filter, "*");
    gtk_file_filter_set_name(filter, "All Files");
    g_list_store_append(filters, filter);

    gtk_file_dialog_set_filters(dialog, G_LIST_MODEL(filters));
    gtk_file_dialog_save(dialog, NULL, NULL, pathFromSaveDialogStringCallback, &destinationPath);
#endif
}

typedef void (*loadingFunction)(std::string);
static void loadDataCallback(GObject *gobject, GAsyncResult *result, gpointer data) {
    GError *error = nullptr;
    GFile *file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(gobject), result, &error);
    if(error == nullptr) {
        std::string path(g_file_get_path(file));

        loadingFunction loadingFunctionPointer = reinterpret_cast<loadingFunction>(data);
        (loadingFunctionPointer)(path);
    }
}

inline void loadFromLoadDialog(loadingFunction loadingFunctionPointer) {
    GtkFileDialog *dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(dialog,
                         NULL,
                         NULL,
                         loadDataCallback,
                         reinterpret_cast<gpointer>(loadingFunctionPointer));
}

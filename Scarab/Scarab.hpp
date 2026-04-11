#include "gtk/gtk.h"
bool update_display();
void handle_run_button();
void draw_spectrum(GtkDrawingArea *area, cairo_t *cr, int width, int height, gpointer data);
void draw_spectrum_frequency(GtkDrawingArea *area,
                             cairo_t *cr,
                             int width,
                             int height,
                             gpointer data);
void draw_spectra_frequency(GtkDrawingArea *area,
                            cairo_t *cr,
                            int width,
                            int height,
                            gpointer data);
void draw_interference_spectrum(GtkDrawingArea *area,
                                cairo_t *cr,
                                int width,
                                int height,
                                gpointer data);
void draw_interference_spectrum_time(GtkDrawingArea *area,
                                     cairo_t *cr,
                                     int width,
                                     int height,
                                     gpointer data);
void draw_interference_phase(GtkDrawingArea *area,
                             cairo_t *cr,
                             int width,
                             int height,
                             gpointer data);
void draw_interference_group_delay(GtkDrawingArea *area,
                                   cairo_t *cr,
                                   int width,
                                   int height,
                                   gpointer data);
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

// Scarab.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "api/OceanDirectAPI.h"
#include "LightwaveExplorerGTK/LightwaveExplorerGraphicalClasses.h"


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

void destroyMainWindowCallback();

bool updateDisplay();
void handleRunButton();
void drawSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void independentPlotQueue();
void handleGetOverlay0();
void handleGetOverlay1();
void handleGetOverlay2();
void handleDeleteOverlay0();
void handleDeleteOverlay1();
void handleDeleteOverlay2();
void handleSave();
void svgCallback();
class OceanSpectrometer {
    std::vector<double> readBuffer;
    std::vector<double> wavelengthsBuffer;
    std::vector<double> overlay0;
    std::vector<double> overlay1;
    std::vector<double> overlay2;

    int lastOverlay = -1;
    int pixelCount = 0;
    bool isInitialized = false;
    long deviceID = 0;
    int error = 0;

public:
    bool hasOverlay0 = false;
    bool hasOverlay1 = false;
    bool hasOverlay2 = false;
	void init(long deviceIDinput) {
		deviceID = deviceIDinput;
		pixelCount = odapi_get_formatted_spectrum_length(deviceID, &error);
        readBuffer = std::vector<double>(pixelCount);
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
    }

    void acquireOverlay(int overlayIndex) {
        std::vector<double> newOverlay(pixelCount);
        odapi_get_formatted_spectrum(deviceID, &error, newOverlay.data(), pixelCount);
        switch (overlayIndex) {
        case 0:
            overlay0 = newOverlay;
            lastOverlay = 0;
            hasOverlay0 = true;
            break;
        case 1:
            overlay1 = newOverlay;
            lastOverlay = 1;
            hasOverlay1 = true;
            break;
        case 2:
            overlay2 = newOverlay;
            lastOverlay = 2;
            hasOverlay2 = true;
            break;
        }
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
            return overlay0.data();
        case 1:
            return overlay1.data();
        case 2:
            return overlay2.data();
        default:
            return nullptr;
        }
    }

    void deleteOverlay(int overlayIndex) {
        switch (overlayIndex) {
        case 0:
            hasOverlay0 = false;
        case 1:
            hasOverlay1 = false;
        case 2:
            hasOverlay2 = false;
        }
    }

    double* data() {
        return readBuffer.data();
    }

    double* wavelengths() {
        return wavelengthsBuffer.data();
    }

    size_t size() {
        return readBuffer.size();
    }

    int getErrorCode() {
        return error;
    }
};
std::vector<OceanSpectrometer> spectrometerSet;

//Main class for controlling the interface
class mainGui {
    bool queueUpdate;
    bool queueSliderUpdate;
    bool queueSliderMove;
    bool queueInterfaceValuesUpdate;
    bool isRunningLive;
    int sliderTarget;
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
    size_t pathTarget;
    int saveSVG = 0;
    bool loadedDefaults = false;
    unsigned int timeoutID;
    mainGui() : queueUpdate(0),
        queueSliderUpdate(0),
        queueSliderMove(0),
        queueInterfaceValuesUpdate(0),
        isRunningLive(false),
        sliderTarget(0),
        pathTarget(0),
        saveSVG(0),
        loadedDefaults(0),
        timeoutID(0) {}
    ~mainGui() {
    }
    void requestPlotUpdate() {
        queueUpdate = true;
    }

    void requestLive() {
        isRunningLive = true;
    }

    void stopLive() {
        isRunningLive = false;
    }

    [[nodiscard]] constexpr bool runningLive() { return isRunningLive; }
    void applyUpdate() {
        if (queueUpdate) {
            queueUpdate = false;
            drawBoxes[0].queueDraw();
        }

        //progressBarBox.queueDraw();

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
        int textWidth = 3;
        int labelWidth = 2;
        int plotWidth = 12;
        int plotHeight = 6;
        int pathChars = 40;
        int colWidth = labelWidth + 2 * textWidth;
        int textCol1a = labelWidth;
        int textCol2a = textCol1a + 2 * textWidth + labelWidth;
        int textCol1b = textCol1a + textWidth;
        int textCol2b = textCol2a + textWidth;
        int buttonCol1 = 3*textWidth;
        int buttonCol2 = buttonCol1 + buttonWidth;
        int buttonCol3 = buttonCol2 + buttonWidth;
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


        textBoxes[0].init(parentHandle, 2*textWidth, 1, textWidth, 1);
        textBoxes[0].setLabel(-2*textWidth, 0, "Integration time (ms)");
        buttons[0].init(("Run"), parentHandle, buttonCol1, 1, buttonWidth, 1, handleRunButton);
        buttons[1].init(("Set O1"), parentHandle, buttonCol1, 3, buttonWidth/2, 1, handleGetOverlay0);
        buttons[2].init(("Set O2"), parentHandle, buttonCol1, 4, buttonWidth / 2, 1, handleGetOverlay1);
        buttons[3].init(("Set O3"), parentHandle, buttonCol1, 5, buttonWidth / 2, 1, handleGetOverlay2);
        buttons[3].init(("Save"), parentHandle, buttonCol1, 6, buttonWidth, 1, handleSave);

        buttons[4].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 3, buttonWidth / 2, 1, handleDeleteOverlay0);
        buttons[5].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 4, buttonWidth / 2, 1, handleDeleteOverlay1);
        buttons[6].init(("\xf0\x9f\x97\x91\xef\xb8\x8f"), parentHandle, buttonCol1+buttonWidth/2, 5, buttonWidth / 2, 1, handleDeleteOverlay2);



        textBoxes[1].init(parentHandle, 0, 2, textWidth, 1);
        textBoxes[2].init(parentHandle, textWidth, 2, textWidth, 1);
        textBoxes[3].init(parentHandle, 2*textWidth, 2, textWidth, 1);

        textBoxes[4].init(parentHandle, 0, 3, textWidth, 1);
        textBoxes[5].init(parentHandle, textWidth, 3, textWidth, 1);
        textBoxes[6].init(parentHandle, 2 * textWidth, 3, textWidth, 1);

        textBoxes[7].init(parentHandle, 0, 4, textWidth, 1);
        textBoxes[8].init(parentHandle, textWidth, 4, textWidth, 1);
        textBoxes[9].init(parentHandle, 2 * textWidth, 4, textWidth, 1);

        textBoxes[10].init(parentHandle, 0, 5, textWidth, 1);
        textBoxes[11].init(parentHandle, textWidth, 5, textWidth, 1);
        textBoxes[12].init(parentHandle, 2 * textWidth, 5, textWidth, 1);

        filePaths[0].init(parentHandle, 0, 6, 4 * textWidth, 1);
        filePaths[0].setMaxCharacters(pathChars);
        filePaths[0].overwritePrint(Sformat("DefaultOutput.txt"));
        
        textBoxes[48].init(window.parentHandle(4), 2, 0, 2, 1);
        textBoxes[49].init(window.parentHandle(4), 4, 0, 2, 1);
        textBoxes[50].init(window.parentHandle(4), 8, 0, 2, 1);
        textBoxes[51].init(window.parentHandle(4), 10, 0, 2, 1);
        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawSpectrum);
        checkBoxes[1].init(("Log"), window.parentHandle(4), 13, 0, 1, 1);
        buttons[10].init(("xlim"), window.parentHandle(4), 0, 0, 1, 1, independentPlotQueue);
        buttons[10].setTooltip("Apply the entered x limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[10].squeeze();
        buttons[11].init(("ylim"), window.parentHandle(4), 6, 0, 1, 1, independentPlotQueue);
        buttons[11].setTooltip("Apply the entered y limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[11].squeeze();
        buttons[12].init(("SVG"), window.parentHandle(3), 5, 0, 1, 1, svgCallback);
        buttons[12].setTooltip("Generate SVG files of the four line plots, with filenames based on the base path set above");
        buttons[12].squeeze();
        console.init(window.parentHandle(1), 0, 0, 1, 1);
        console.cPrint("Attached spectrometers:\n");
        
        g_signal_connect(window.window, "destroy", G_CALLBACK(destroyMainWindowCallback), NULL);
        window.present();
        initializeSpectrometers();
        pulldowns[0].init(parentHandle, buttonCol1, 0, buttonWidth, 1);
        pulldowns[0].setLabel(-labelWidth, 0, "Device:");
        drawBoxes[0].queueDraw();
        timeoutID = g_timeout_add(50, G_SOURCE_FUNC(updateDisplay), NULL);



        //getSpectrum();
    }

    void initializeSpectrometers() {
        int deviceCount = odapi_probe_devices();
        int deviceIdCount = odapi_get_number_of_device_ids();
        std::vector<long> deviceIds(deviceIdCount);
        int retrievedIdCount = odapi_get_device_ids(deviceIds.data(), deviceIdCount);
        int error = 0;
        spectrometerSet = std::vector<OceanSpectrometer>(deviceCount);
		for (int i = 0; i < deviceCount; i++) {
			odapi_open_device(deviceIds[0], &error);
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


void independentPlotQueue() {
    //theGui.requestPlotUpdate();
    //theGui.applyUpdate();
    theGui.console.cPrint("works\n");
}

void handleGetOverlay0() {
    spectrometerSet[0].acquireOverlay(0);
}

void handleGetOverlay1() {
    spectrometerSet[0].acquireOverlay(1);
}

void handleGetOverlay2() {
    spectrometerSet[0].acquireOverlay(2);
}

void handleDeleteOverlay0() {
    spectrometerSet[0].deleteOverlay(0);
}
void handleDeleteOverlay1() {
    spectrometerSet[0].deleteOverlay(1);
}
void handleDeleteOverlay2() {
    spectrometerSet[0].deleteOverlay(2);
}

void handleSave() {

}


void drawSpectrum(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
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
    
    if (theGui.runningLive()) {
        spectrometerSet[activeSpectrometer].setIntegrationTime(1000 * theGui.textBoxes[0].valueDouble());
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
    sPlot.dx = 1;
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
    sPlot.forcedXmin = yMax;
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
        if (sPlot.ExtraLines > 1 && firstAdded != 1 && spectrometerSet[activeSpectrometer].hasOverlay1) sPlot.data3 = spectrometerSet[0].getOverlay(1);
        else if (sPlot.ExtraLines > 1 && firstAdded != 2 && spectrometerSet[activeSpectrometer].hasOverlay2) sPlot.data3 = spectrometerSet[0].getOverlay(2);
       
        if (sPlot.ExtraLines > 2) sPlot.data4 = spectrometerSet[0].getOverlay(2);
    }
    sPlot.plot(cr);

    
}
bool updateDisplay() {
    theGui.requestPlotUpdate();
    theGui.console.updateFromBuffer();
    theGui.applyUpdate();
    return true;
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
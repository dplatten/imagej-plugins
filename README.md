This repository contains my ImageJ plugins to help with some medical physics tasks. The skeleton I have used was taken from https://github.com/imagingbook/imagej1-plugins-ide-setup. I have stripped out the netbeans and Eclipse folders.

# The plugins

## TextureMap
Diagnostic Radiology. Calculate a variance or signal-to-noise map of a
digital x-ray image.

The routine will run on a user-selected rectangular region of
interest (ROI). If no ROI is selected then the whole image is
selected and analysed.

It breaks the large ROI into smaller square ROIs, and calculates
the chosen statistical parameter for each region. The user is
asked to supply the side-length of the smaller ROIs (in mm).

The routine calculates variance as:

sum[x - mean(x)] / n

The routine calculates standard deviation as:

sqrt(variance)

The routine calculates signal-to-noise as:

mean(x) / standard deviation(x)

## ResampleImage
Nuclear Medicine. This re-samples an existing image to simulate a new image that would
have been acquired using a lower number of counts.

## TaskTransferFunction
Diagnostic Radiology. This calculates the task transfer function of a circular
object within an image.

## DetectabilityIndex
Uses the TaskTransferFunction, NoisePowerSpectrum, human eye model and task function
to calculate a numerical index that relates to the detectability of an object. Note
that only the TaskTransferFunction part of this is written so far.


# Using this repository

## Setup for writing ImageJ Plugins with IntelliJ
This repo contains a minimal setup for writing ImageJ (1) plugins with [IntelliJ](https://www.jetbrains.com/idea/).
The project is set up with ``<project-root>/plugins/`` as the default output folder (for generated ``.class`` files).

This repository is based on the [**imagingbook**](http://imagingbook.com) support suite.
See [www.imagingbook.com](http://imagingbook.com) for additional resources.

### Setup
Clone this repository. It contains a folder for ``project-intellij/``, a self-contained project for IntelliJ.

#### IntelliJ:
Start the **IntelliJ IDEA** and use ``Open`` in the *Welcome screen* to navigate to the ``project-intellij/`` folder.
Editing and saving plugin source files should update the associated class files (in ``plugins/``) automatically.

### Starting ImageJ
The ImageJ runtime can be launched in various ways:
- **Windows**: Execute ``ImageJ.exe`` (by double-clicking on the file).
When ImageJ starts up, it may ask for the ``javaw.exe`` executable, typically located in ``C:\Program Files\java\jre1.8xxx\bin\``. In case of problems, simply delete the ``ImageJ.cfg`` file and start anew.
- **MacOS**: Launch ``ij.jar``.
- **Java**: Run the ``ij.ImageJ.main()`` method within Eclipse.

The entire ImageJ functionality is contained in the single archive ``ij.jar``. To **update** to the most recent version, simply select ``Help`` -> ``Update ImageJ...`` from the ImageJ main menu.

### Adding/editing your plugin code
Code for ImageJ plugins is contained in the ``<project-root>/src-plugins/`` directory. Plugins may be contained in Java packages (such as ``my_plugins`` in this example). **Note** that packages with plugins may only be **one level deep**, otherwise ImageJ will not find them! It is recommended to use at least one underscore (``_``) in a plugin name to make ImageJ automatically install the plugin into the ``Plugins`` menu at startup.

### Executing plugins
At startup, ImageJ automatically installs existing plugins (under the above conditions) into the ``Plugins`` menu. To execute, simply select the listed plugin from the menu.

When the plugin's source code is **edited** in the IDE, the associated ``.class`` file in ``plugins/`` is updated (by Eclipse/IntelliJ), but **not** automatically reloaded by the ImageJ runtime. To **exectute an edited plugin** in ImageJ, use ``Plugins`` -> ``Compile and Run...`` and select the associated ``.class`` file (no compiler is needed).

### Adding other libraries (jars)
This project uses **no dependency management** (such as *Maven*) to keep things simple. If any external libraries are required, just do the following:
- Copy the associated JAR file ``xxx.jar`` into ``jars/``.
- In your IDE, add the JAR file to Java build path.
- Restart ImageJ.


### Additional ImageJ resources

- [ImageJ Home](https://imagej.nih.gov/ij/index.html)
- [ImageJ Plugins](http://rsbweb.nih.gov/ij/plugins/index.html)
- [ImageJ API (JavaDoc)](http://rsbweb.nih.gov/ij/developer/api/index.html)
- [ImagingBook (book and source code)](http://imagingbook.com)

// @(#)root/ged:$Id$
// Author: Denis Favre-Miville   08/09/05

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "HelpSMText.h"

const char gHelpSMTopLevel[] =
"The Style Manager handles a list of styles in the ROOT session. It loads\n"
"by default the five styles provided by ROOT: Default, Plain, Bold, Video,\n"
"and Pub. If this list does not contain a style with the characteristics\n"
"you want, you can create a new one and than apply it.\n"
"\n"
"The Style Manager interface is composed of two parts:\n"
"   - the top level interface that manages a list of styles;\n"
"   - the style editor, which deals with the current style settings.\n"
"\n"
"'Available Styles':\n"
"   Contains the list of available styles for the current ROOT session and the\n"
"   currently selected one. The field on the right shows the setting of the\n"
"   gStyle. You can set the global variable gStyle to the selected style by the\n"
"   button in the middle.\n"
"\n"
"'Apply on':\n"
"   Displays information for the currently selected canvas and object in the\n"
"   ROOT session. This selection might be changed by clicking on another\n"
"   object with the middle mouse button. You have a choice to apply the\n"
"   selected style on the selected object or on all available canvases.\n"
"   in the ROOT session.\n"
"   WARNING: You cannot undo the changes after applying the style! If you are\n"
"   not sure of that action, it may be better to see a preview of what you are.\n"
"   going to apply.\n"
"\n"
"Preview\n"
"   Shows a preview of the selected canvas according to the selected style.\n"
"\n"
"Run Time Preview\n"
"   Updates the preview any time a value of the selected style is changed.\n"
"   For drawings that take a time it is better to disable this option.\n"
"\n"
"Create a new style / delete a style:\n"
"   You can access these functionalities via the menu or the tool bar.\n"
"   If you create a style, a clone of the selected style will be created;\n"
"   you will have the opportunity to modify it later via the editor.\n"
"   Moreover, during the creation process, you will be asked for a name and\n"
"   a description for your brand new style. The name can obviously not be the\n"
"   same as another already existing style.\n"
"   When you choose 'delete', the selected style is removed from the ROOT\n"
"   session. The selected style will be lost if you didn't saved it in a C++\n"
"   macro file before. Nota Bene: You are not allowed to delete gStyle.\n"
"\n"
"Export a style (in a C++ macro file) / Import a style (from a macro):\n"
"   To store a style longer than for one session of ROOT or in order to\n"
"   share some styles with others, you can save it in a C++ macro file.\n"
"   This can be done via the menu and the tool bar. The selected style will\n"
"   be saved. The name of the macro must be 'Style_*.C', where * can be\n"
"   replaced by anything you want.\n"
"   A style's macro can be imported at any time. The new style created\n"
"   in ROOT will become the selected style. WARNING: When a style is loaded\n"
"   from a macro, every already existing style with the same name will be\n"
"   overwritten.\n"
"\n"
"Apply a style (on an object) / Import a style (from an object):\n"
"   There is a specific button to apply a style on an object and its use\n"
"   has already been described in the overview.\n"
"   To import a style from a canvas, the 'Style / Import from... / Canvas'\n"
"   menu is the only way. A new style will be created in the ROOT session and\n"
"   will become the selected style. This style is a clone of gStyle where\n"
"   every style's information contained in the selected canvas (the canvas\n"
"   containing the current selected pad) is written. You can consequently\n"
"   import a style from a canvas and apply it later on another object.\n"
"\n"
"Preview a style:\n"
"   The predicted result when applying a style can be seen if the 'Preview'\n"
"   checkbutton is selected. The preview includes the original canvas.\n"
"\n"
"Editor's buttons:\n"
"   Open / close the editor:\n"
"      A specific button to open or close the editor exists in the style\n"
"      manager. Its caption is 'Edit >>' when the editor is closed and\n"
"      switches to 'Close <<' when it is opened.\n"
"   Reset a style (to a previously saved state):\n"
"      When the editor is opened, in the bottom of the main window of the\n"
"      style manager, a 'Reset' button enables you to reset the values of the\n"
"      selected style. So doing, you cancel all changes made since the last\n"
"      time you saved that style in a macro file. If the selected style is\n"
"      one of the five predefined styles of ROOT (Plain, Bold, Video, Pub or\n"
"      Default), it will be reset using the specific code in the ROOT files.\n"
"   Update the preview:\n"
"      An 'Update Preview' is also available when a preview is shown and\n"
"      the run time option is disabled. This button allows you to refresh the\n"
"      preview at any time.\n"
;

const char gHelpSMGeneral[] =
"Sets the general style values used when creating new objects\n"
"\n"
"'Fill' group: \n"
"   Sets the default filling color and pattern for drawing of any new object;\n"
"   the hatchings' line width and spacing for patterns > 3100.\n"
"\n"
"'Text' group:\n"
"   Sets the text color, font, alignment, font size and angle.\n"
"   The font size can be specified in pixels or in % of the pad size.\n"
"\n"
"'Line' group:\n"
"   Sets the parameters for creating a new line: color, width and line style.\n"
"\n"
"'Marker' group:\n"
"   Sets the default parameters markers: color, style and size. \n"
"\n"
"Screen factor:\n"
"   Sets a coefficient for different screen resolutions.\n"
"   For exemple, if a macro creates a 800*600 canvas, and the screen\n"
"   factor is set to 0.8, a 640*480 canvas will be shown.\n"
;

const char gHelpSMCanvas[] =
"Sets the style values related to the canvas objcts.\n"
"\n"
"'Fill' group:\n"
"   Sets the canvas fill color.\n"
"\n"
"'Geometry' group:\n"
"   Sets the canvas geometry:the width, the height, and the (x,y) position \n"
"   where any new canvas will be placed, by default.\n"
"\n"
"'Border' group:\n"
"   Sets the canvas border mode as sunken, none or raised; and the border\n"
"   width. These parameters are closely linked to the canvas fill color\n"
"   because of the 3D effect will be based on it.\n"
"\n"
"'Date' group:\n"
"   You can choose to show or hide the date in your canvases via the\n"
"   check button 'Show'.\n"
"   If 'Show' is selected, you can set the following date parameters:\n"
"   color, size, format, font, alignement and angle. The date position \n"
"   is specified in percent of the pad size.\n"
;

const char gHelpSMPad[] =
"Sets the style values related to the pad objcts.\n"
"\n"
"'Margin' group:\n"
"   Generally, when including an object in a pad, one doesn't want it to\n"
"   fill the whole pad. These margins describe the amount of empty space left\n"
"   around the objects inside a pad.\n"
"\n"
"'Border' group:\n"
"   Sets the pad border mode as sunken, none or raised; and the border\n"
"   width. These parameters are closely linked to the pad fill color\n"
"   because of the 3D effect will be based on it.\n"
"\n"
"'Fill' group: \n"
"   You can modify here the color which is used to fill the pads.\n"
"\n"
"'Ticks' group:\n"
"   This group allows you to show the special pad ticks along each axis.\n"
"\n"
"'Grid' group:\n"
"   You can choose to show or hide the grid along the two axis X and Y, by\n"
"   selecting the options. Moreover, you can change the lines settings:\n"
"   the color, the line width and line style can be modified.\n"
;

const char gHelpSMHistos[] =
"Sets the style values related to the histogram objcts. The style values\n"
"are grouped in three sub-tabs: Histos, Frames and Graphs.\n"
"\n"
"'Histos' Tab:\n"
"\n"
"'Fill' group: \n"
"   Sets the color and the pattern for filling histograms.\n"
"\n"
"'Line' group:\n"
"   Sets the line parameters for histograms: color, line width and style\n"
"\n"
"'Bar' group:\n"
"   Sets the width and offset of bars for graph objects.\n"
"\n"
"'Contours' group:\n"
"   Sets the number of contours for drawing by default.\n"
"\n"
"'Axis' group:\n"
"    'Minimum zero':\n"
"    Check this option to force the minimum of axis to be zero.\n"
"    Otherwise, the minimum will be computed.\n"
"    'Paint format':\n"
"    It is the axis printing format for the TH2 objects. This format is\n"
"    used as an argument in a C++ printf() method.\n"
"\n"
"'3D Cylindrical' group:\n"
"   The inner radius is useful for some 3D graphics. Using the cylindrical\n"
"   coordinates system, this value represents the percentage of the radius\n"
"   allocated to the inner radius.\n"
"\n"
"The 'Frames' Tab:\n"
"'Fill' group: \n"
"   Sets the color and the pattern for frames.\n"
"\n"
"'Line' group:\n"
"   Sets the line parameters for frames: color, line width and style\n"
"\n"
"'Border' group:\n"
"   Sets the frame border mode as sunken, none or raised; and the border\n"
"   width. These parameters are closely linked to the frame fill color\n"
"   because of the 3D effect will be based on it.\n"
"\n"
"\n"
"The 'Graphs' tab:\n"
"'Line' group:\n"
"   Sets the line parameters for graphs: color, line width and style\n"
"\n"
"'Errors' group:\n"
"    End error size:\n"
"      Specifies the size of lines which are drawn at the end of error\n"
"      bars.\n"
"    Error along X:\n"
"      Sets percent of the bin width to be used for errors along X.\n"
"\n"
"'Draw Border' option:\n"
"   Specifies if the border of the filled functions must be drawn or not.\n"
;

const char gHelpSMAxis[] =
"Sets the style values related to the axis grouped in three sub-tabs:\n"
"X axis, Y axis and Z axis.\n"
"\n"
"'Decimal labels' option:\n"
"   Allows you to specify if you want to print the decimal part of the\n"
"   labels in the axis.\n"
"\n"
"'Date/Time offset' group:\n"
"   Sets the time offset shown in two parts: as a date and time.\n"
"\n"
"'Line' group:\n"
"   Sets the color, the ticks' length and logarithmic scale.\n"
"\n"
"'Title' group:\n"
"   Sets the appearance of the axis title: color, font, font size and offset.\n"
"\n"
"'Divisions' group:\n"
"    * Sets the tertiary, the secondary and the primary axis division numbers.\n"
"    * with the 'Optimize' option ON or OFF\n"
"\n"
"'Labels' group:\n"
"   Sets the axis labels' appearance. Similar as the title group.\n"
;

const char gHelpSMTitle[] =
"In this tad, the information related with titles are edited. You can choose\n"
"to show or hide the title in your canvases via the check button 'Show title'.\n"
"\n"
"'Fill' group: \n"
"   You can modify here the color and the pattern which is used to fill\n"
"   the title pave.\n"
"\n"
"'Text' group:\n"
"   The various parameters in the text group are the color, the font, the\n"
"   alignement and the font size. The font size can be specified in pixels or\n"
"   in percentage of the pad.\n"
"\n"
"'Shadow' group:\n"
"   Shadows are little borders added behind the object to create an\n"
"   illusion of third dimension. The size of the shadow depend on the\n"
"   importance of each object.\n"
"\n"
"'Geometry' group:\n"
"   What can be changed here is the geometry of the title pave. By setting\n"
"   these values, you define the size and the location where the title pave\n"
"   will be placed. The fields are the abscissa (X), the ordinate (Y),\n"
"   the width (W) and the height (H).\n"
;

const char gHelpSMStats[] =
"In this tad, the information related with stats are edited.\n"
"\n"
"'Fill' group: \n"
"   You can modify here the color and the pattern which is used to fill\n"
"   the statistics pave.\n"
"\n"
"stats' shadow:\n"
"   Shadows are little borders added behind the object to create an\n"
"   illusion of third dimension. The size of the shadow depend on the\n"
"   importance of the object.\n"
"\n"
"'Text' group:\n"
"   The various parameters in the text group are the color, the font and\n"
"   the font size. The font size can be specified in pixels or in percentage\n"
"   of the pad.\n"
"\n"
"'Geometry' group:\n"
"   What can be changed here is the geometry of the stats pave. By setting\n"
"   these values, you define the size and the location where the stats pave\n"
"   will be placed. The fields are the abscissa (X), the ordinate (Y),\n"
"   the width (W) and the height (H).\n"
"\n"
"'Stat Options' group:\n"
"   You can select the information you want to see in the stats pave.\n"
"   The paint format is used as an argument in a C++ printf() method: it\n"
"   describes how the value will be shown.\n"
"\n"
"'Fit Options' group:\n"
"   You can select the way you want to show information in the pave.\n"
"   The paint format is used as an argument in a C++ printf() method: it\n"
"   describes how the value will be shown.\n"
;

const char gHelpSMPSPDF[] =
"In this tad, the information related with PS and PDF are edited. These\n"
"settings can't be previewed, but as in the other tabs, any modification in\n"
"this one will take effect.\n"
"\n"
"'Header' group:\n"
"   Define a string to be inserted in the Postscript header. The string in\n"
"   header will be added to the Postscript file immediatly following the\n"
"   %%Page line.\n"
"   For example, this string may contain special Postscript instructions\n"
"   like: '200 200 translate'\n"
"   The following header string will print the string \"my annotation\" at\n"
"   the bottom left corner of the page (outside the user area):\n"
"   'gsave 100 -100 t 0 r 0 0 m /Helvetica-Bold findfont 56 sf 0 0 m ( my\n"
"    annotation ) show gr'\n"
"\n"
"'Title' group:\n"
"   Define a string to be used in the %%Title of the Postscript files.\n"
"   If this string is not defined, ROOT will use the canvas title.\n"
"\n"
"'Paper Size' group:\n"
"   Set paper size for PostScript output. You can choose a predefined size\n"
"   or specify manually the width and the height of the printed area.\n"
"\n"
"'Line scale' option:\n"
"   To specify a line scale factor when drawing lines on Postscript.\n"
"\n"
"'Color Model' group:\n"
"   Define the color model used by TPostScript and TPDF (RGB or CMYK).\n"
"   CMYK model is a subtractive color model unlike RGB which is an additive.\n"
"   They are mainly used for printing purposes. CMYK means Cyan Magenta\n"
"   Yellow Black, RGB means Red Green Blue.\n"
;

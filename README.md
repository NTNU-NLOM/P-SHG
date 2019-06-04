# P-SHG  

# # The original code was made by NTNU-NLOM, here : https://github.com/NTNU-NLOM/P-SHG  !!!

This is a FORKED version of the original, to be controlled with the following instruments:  
- a Thorlabs PM100 power-meter OR a National instrument DAQ card
- a Micos Pollux motor 
- a Newport motor controlled by an ESP-100 or 300 controller
- a PRM1/Z8 from Thorlabs

The designed has been slightly changed to allow to fit the screen, and the code modified for the easy use of other motors.
For the NI card, user should be sure to have installed the NI libraries, put the correct paths in the file and installed MinGW64+ (for NI x32 version) as prompted by Matlab.  

UPDATE : - NI card as effect-meter
- possibilities to combine results
- possible to continue calib in the end
- skip some angles

Made by : Max PINSARD, INRS-EMT, Varennes, QC, Canada.  
August 2018 - updated May 2019
maxime.pinsard@outlook.com  

------

The matlab code to run an automated polarization-resolved second harmonic generation (P-SHG) setup designed for a commercial microscope.
For more details on application and implementation see "Romijn, Elisabeth Inge; Finnøy, Andreas; Kumar, Rajesh; Lilledahl, Magnus Borstad. (2018) Automated calibration and control for polarization-resolved second harmonic generation on commercial microscopes. PLoS ONE. vol. 13 (4)".

Two graphical user interfaces (GUIs) were created; one to run the calibration of the P-SHG setup ("polarization_calibration"), and the other to enable P-SHG imaging ("polarization_imaging"). (The GUI "properties_polarization_calibration" is called upon by "polarization_calibration".)


This code is meant as an example, and has to be adapted to accomodate the employed hardware. 

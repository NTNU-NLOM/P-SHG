# P-SHG
The matlab code to run an automated polarization-resolved second harmonic generation (P-SHG) setup designed for a commercial microscope.
For more details on application and implementation see "Automated calibration and control for polarization resolved second harmonic generation on commercial microscopes" (in preparation).

Two graphical user interfaces (GUIs) were created; one to run the calibration of the P-SHG setup ("polarization_calibration"), and the other to enable P-SHG imaging ("polarization_imaging"). (The GUI "properties_polarization_calibration" is called upon by "polarization_calibration".)

The hardware this software is designed for consists of:
- Calibration: Three motorized rotation stages (PRM1/MZ8, Thorlabs) and a power meter (Labmax-TOP, Coherent)
- Imaging: Two motorized rotation stages (PRM1/MZ8, Thorlabs) and a confocal microscope (TCS SP8, Leica)

This code is ment as an example, and has to be adapted to accomodate the employed hardware. 

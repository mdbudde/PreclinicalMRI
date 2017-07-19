# PreclinicalNeuroMRI
A collection of tools for preclinical brain and spinal cord MRI

NOTE: I originally wanted to share source code, which facilitates installation.  However, I realize that this could violate Bruker's license agreement.  I intend to ask Bruker for clarification regard sharing of sequences.  I may need to share the sequences in binary format, but that may complicate the installation since some files get installed during compilation.  I need to followup on these details, but in the meantime, I have not made the source code publically available.


Modified sequences:

gatedMGEv1: 
Uses external respiratory/cardiac gating but maintains steady state (dummy scans) between triggers.  Typically used for 3D multiple-gradient echo images of the spinal cord.

dtiEpi_GRFat: 
Includes a gradient reversal fat saturation (flipping the sign of the refocussing pulse) in conjunction with fat suppression pulse.
Also includes an option for reversed 'blipping' of the EPI readout.

vfaRARE: 
Implements a variable flip angle Fast Spin Echo (Bruker: RARE).  Used for high resolution 3D spin echo imaging.  Allows a rare factor of up to 76 while keeping effective TEs similar to typical RARE.  Can be used with linear or centric encoding.

shapedDwi_Epi: 
A sequence to enable any arbitrary gradient waveform to be implemented for diffusion weighting using an external file.  Compatible with oscillating gradient spin echo (OGSE), rotation of the q vector (qMAS), and many others.
This sequence will require the user to create/modify their own waveform file.  The examples included in the sequence, ShapedDwiStick.gp, etc should be modified to accomodate the users gradient characteristics.  Re-compiling the sequence will install them to the correct directory. Clicking on the 'Reload Shapes' button in the PV interface will ensure they get reloaded, obviously.



Simulations:
The simulation folder contains matlab routines to create "beaded" geometries for use in camino diffusion simulations.


Matthew Budde, PhD
Medical College of Wisconsin

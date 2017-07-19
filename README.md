# PreclinicalMRI
A collection of tools for preclinical brain and spinal cord MRI

NOTE: I originally wanted to share source code, which facilitates installation.  However, I realize that this could violate Bruker's license agreement.  I intend to ask Bruker for clarification regard sharing with other users.  The sequencesi may need to be shared in binary format, but that may complicate the installation.  In the meantime, the modified source code sequences are not publically available.

Please contact me directly with questions or availability of sequences: mdbudde@users.noreply.github.com
 

List of available modified sequences:
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

fFOV_DwiEpi: (Not yet released)
A reduced field of view diffusion EPI sequence using 2D echo planar excitation (2DRF).

dde_press: (Not yet release)
A press spectroscopy sequence with double diffusion encoding preparation.


Simulations:
The simulation folder contains matlab routines to create "beaded" geometries for use in camino diffusion simulations.


Matthew Budde, PhD
Medical College of Wisconsin

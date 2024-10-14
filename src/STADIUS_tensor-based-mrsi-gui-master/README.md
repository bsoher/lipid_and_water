**Tensor-based MRSI GUI**
**Introduction**
1.	Application
MRSI signal processing.
2.	Purpose
Tensor-based approach to residual water suppression and tissue type differentiation using Magnetic Resonance Spectroscopic Imaging (MRSI) signals. 
3.	Description
This MATLAB based graphical user interface (GUI) implements tensor-based approaches for residual water suppression and tissue type differentiation using Magnetic Resonance Spectroscopic Imaging (MRSI) signals. The algorithms for residual water suppression and tissue type differentiation are described in [1] and [2], respectively.
4.	Citations
[1] Bharath, H. N., Debals, O., Sima, D. M., Himmelreich, U., De Lathauwer, L., and Van Huffel, S. Tensor based method for residual water suppression in 1H magnetic resonance spectroscopic imaging. IEEE Transactions on Biomedical Engineering, vol. 66, no. 2 (Feb. 2019), pp. 584-594. DOI: 10.1109/TBME.2018.2850911
[2] Bharath, H. N., Sima, D. M., Sauwen, N., Himmelreich, U., De Lathauwer, L., and Van Huffel, S. Nonnegative canonical polyadic decomposition for tissue-type differentiation in gliomas. IEEE Journal of Biomedical and Health Informatics 21, 4 (July 2017), 1124–1132.
5.	Keywords/tags
Magnetic resonance spectroscopic imaging, Löwner matrix, Hankel matrix, Blind source separation, tensor decomposition

**Describing The Code**
Contents: (Main Functions)
1.	MRSI_tissue.m: Matlab code for the GUI
2.	MRSI_tissue.fig : the figure for the GUI, no need to open this one
3.	MRSI_GUI_user_manual.pdf: the documentation
4.	CSI_spectro_35_SENSE_11_2_raw_act. SDAT/SPAR: Example MRSI data from Philips scanner.
5.	aa_rep_img.JPG: Example background MRI image. 
6.	Loewner_water_supression.mat: Function to suppress residual water using Löwner-tensor based method.
7.	water_removal_HSVD.mat: Function to suppress residual water using HSLVD method.
8.	Hankel_water_supression.m: Function to suppress residual water using hankel-tensor based method.
The code (MRSI_tissue.m) implements a graphical user interface for visualizing MRSI signals, residual water suppression in MRSI and performing tissue type differentiation. The code ‘Loewner_water_supression.m’, ‘water_removal_HSVD.m’ and ‘Hankel_water_supression’ can be used to suppress residual water in MRSI without the GUI. 

**Running The Code**
To run the GUI, simply type MRSI_tissue in the Matlab command line and press enter. The GUI will initialize.
Basic GUI operation pipeline: 
1)	load the MRSI data (eg: CSI_spectro_35_SENSE_11_2_raw_act. SDAT/SPAR)
2)	Load the background MRI image (eg: aa_rep_img.JPG). GUI will also work without loading the background MRI image.
3)	Crop the voxels outside the volume of interest using crop button.
4)	Perform residual water suppression using HLSVD/Hankel/ Löwner button.
5)	Perform tissue differentiation using NCPD button.

**Program Output**
The outputs of the program are displayed in the GUI window. It can be saved in matlab data format and loaded again by the GUI to view the outputs. 

**Others**
The software expects that background MRI image slice is co-registered and aligned with the MRSI voxel grid. 

**License**
See the Software License Agreement file for license rights and limitations. By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and condition as specified in the License agreement.


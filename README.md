# 3D-Divergence-Free-Flow-Reconstruction
This folder contains Matlab code of 3D flow divergence free reconstruction using 3D distribution of 3D, 2D or even 1D velocity vector components input. The input velocity vectors can be acquired from ultrasound Doppler (1D to 4D), Ultrasound speckle tracking (usually 2D to 4D), MRI Phase Contrast imaging (1D to 4D), Particle Imaging Velocimetry (usually 2D to 4D), etc.
It assumes flows are incompressible, and to construct 3D full field flow it requies scattered 3D, 2D or 1D velocity components to spatially fill the 3D interested volume:
1. with 3D distribution of 3D velocity component input. It is usually a well-posed problem and is symmemtric, and can be solved by e.g., Conjugate gradient method.
2. with 3D distribution of 2D velocity component input. It assumes 2D velocities are planar projection of 3D flow velocities onto the imaging plane, and conduct inverse with regularization. The author uses GMRES with early-stopping to solve the system.
3. with 3D distribution of 1D velocity component input. It assumes 1D velocities are projection of 3D flow velocities onto the imaging beam direction (which is a vector), and conduct inverse with regularization. The author uses GMRES with early-stopping to solve the system.

Hyperparameters to tune for efficiency and accuracy:
1. Shape parameter
2. Regularization parameter, i.e., the early stopping convergence criteria

To use this code, you must cite:
1) Zhou X, Papadopoulou V, Leow CH, Vincent P, Tang M-Xet al., 2019, 3-D flow reconstruction using divergence-free interpolation of multiple 2-D contrast-enhanced ultrasound particle imaging velocimetry measurements, Ultrasound in Medicine and Biology, Vol: 45, Pages: 795-810, ISSN: 0301-5629
2) Zhou X, Vincent P, Zhou X, Leow CH, Tang M-Xet al., 2019, Optimization of 3-D Divergence-Free Flow Field Reconstruction Using 2-D Ultrasound Vector Flow Imaging, Ultrasound in Medicine and Biology, Vol: 45, Pages: 3042-3055, ISSN: 0301-5629
3) Zhou X, Zhou X, Leow CH, Vincent P, Tang Met al., 2018, 3D Flow Reconstruction and Wall Shear Stress Evaluation with 2D Ultrafast Ultrasound Particle Imaging Velocimetry, IEEE International Ultrasonics Symposium (IUS), Publisher: IEEE, ISSN: 1948-5719

With any queries, contact the author: zhouxinhuan0205@126.com.

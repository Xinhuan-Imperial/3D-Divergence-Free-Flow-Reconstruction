# 3D-Divergence-Free-Flow-Reconstruction
This folder contains Matlab code of 3D flow divergence free reconstruction using 3D distribution of 3D or 2D velocity vectors input. The input velocity vectors can be acquired from ultrasound Doppler, Ultrasound speckle tracking, MRI Phase Contrast imaging, Particle Imaging Velocimetry, etc.
It assumes flows are incompressible, and to construct 3D full field flow it requies 3D or 2D velocites to fill the 3D interested volume:
1. with 3D velocity input. It is usually a well-posed problem and is symmemtric, and can be solved by e.g., Conjugate gradient method.
2. with 2D velocities input. It assumes 2D velocities are planar projection of 3D flow velocities onto the imaging plane, and conduct inverse with regularization. The author uses GMRES with early-stopping to solve the system.

To use this code, you must cite:
1) Xinhuan Zhou, Virginie Papadopoulou, Chee Hau Leow, Peter Vincent, Mengxing Tang. 3D Flow Reconstruction Using Divergence Free Interpolation of Multiple 2D Contrast Enhanced Ultrasound Particle Imaging Velocimetry Measurements. Ultrasound in Medicine & Biology
2) Xinhuan Zhou, Peter Vincent, Xiaowei Zhou, Chee Hau Leow, Mengxing Tang. Optimization of a Fast 3D Divergence Free Flow Reconstruction Using 2D Ultrasound Vector Flow Imaging. Ultrasound in Medicine & Biology
With any queries, contact the author: zhouxinhuan0205@126.com.

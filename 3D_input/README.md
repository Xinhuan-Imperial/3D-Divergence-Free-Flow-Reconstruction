Input: 3D sparse/scattered distribution of 3D velocity component vectors. The interested 3D flow domain has to be filled by such 3D sparse velocity vectors. The best accuracy is achieved with uniform as well as dense 3D distribution of velocity vectors. Only Matlab code is given.

Output: 3D interpolated velocity vectors (can be full field flow in the domain), with the assumption that flow being incompressible.

Parameters to tune: 

1. shape parameter kernel.alpha (see the example code), usually set to be reciprocal of the average Eucledian distance between two nearest scattered input velocity vectors. Small epsilon--> instability to input error; Large epsilon--> low fitting accuracy. Uniformly distributed velocity vectors in space is preferred, and epsilon can be 1/aver_dist. However usually an intermediate shape parameter results in highest accuracy.

2. regularization parameter tol (see the example code), usually set to be the relative error of input velocity vectors, as the input may contain measurement/registration error etc. However usually an intermediate regularization parameter results in highest accuracy.

Examples:
1. Poiseulle_flow_3dinput.m: steady Poiseulle flow in a staight pipe, diameter 5mm, length 10mm. Input is from analytical solutions.
2. Womersley_flow_3dinput.m: unsteady Womersley flow in a staight pipe, diameter 5mm, length 10mm.  Input is from analytical solutions.

To try with velocity input from real measurements that contain high level of noise, you must tune shape parameter and regularization parameters very carefully for high accuracy.

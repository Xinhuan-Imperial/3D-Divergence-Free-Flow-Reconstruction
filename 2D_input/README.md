Input: 2D sparse/scattered velocity vectors. The interested 3D flow domain has to be filled by such 2D projected velocity vectors.

Output: 3D interpolated velocity vectors (can be full field flow in the domain), with the assumption that: 1. flow being incompressible; 2. 2D velocity vectors are planar projection of 3D velocity vectors onto the imaging plane.

Software: Matlab

Parameters to tune: 
1. shape parameter epsilon, usually set to be average distance between two nearest scattered input velocity vectors
2. stopping criteria, usually slightly larger than the relative error of input 2D velocity vectors. For example, if the input velocities are acquired by Ultrasound speckle tracking, with a relative error of 10% (consider all error sources), the stopping criteria can be set as 1.1*10%, i.e., 0.11. You may need to tune different parameters for a balance between fitting accuracy and error robustness.

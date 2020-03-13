Input: 3D sparse/scattered velocity vectors. The interested 3D flow domain has to be filled by such 3D sparse velocity vectors. The best accuracy is achieved with uniform 3D distribution of velocity vectors.

Output: 3D interpolated velocity vectors (can be full field flow in the domain), with the assumption that flow being incompressible.

Software: Matlab

Parameters to tune: shape parameter epsilon, usually set to be reciprocal of the average Eucledian distance between two nearest scattered input velocity vectors. Small epsilon--> instability to input error; Large epsilon--> low fitting accuracy. Uniformly distributed velocity vectors in space is preferred, and epsilon can be 1/aver_dist.

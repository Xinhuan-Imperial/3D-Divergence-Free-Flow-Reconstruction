%% by Dr. Xinhuan Zhou from Imperial College London, zhouxinhuan0205@126.com
% This code is for reconstruction of steady Poiseulle flow in a ciruclar pipe
clear;close all;clc;

%% parameters of the straight pipe and input velocity vector distribution
cylinder.length_limit=[0 10];%total length along z-axis
cylinder.center=[0 0];
cylinder.radius=5;
Q=10000;%volumetric flow rate

%spacing along two in-plane and out-of-plane translation directions
x_w=1;
y_w=1;
spacing=1;
scale=4;%scaling parameter for Matlab quiver function

%%
points=[];
[x,y]=meshgrid(-cylinder.radius:x_w:cylinder.radius,-cylinder.radius:y_w:cylinder.radius);
idx=x.^2+y.^2<=cylinder.radius^2;
x=x(idx);
y=y(idx);
for z=cylinder.length_limit(1):spacing:cylinder.length_limit(2)
    points=[points;x(:),y(:),z*ones(numel(x),1)];
end
vel_in=v_poiseulle(points,cylinder,Q);
colored_quiver(points(:,1),points(:,2),points(:,3),vel_in(:,1),vel_in(:,2),vel_in(:,3),scale,' ');
title('input 3D velocity');

%% construct the DFI matrix
stopping_constant=1.05;%early-stopping constant
tol=1e-3;%input is free of error, so tol can be very small here. If error exists, increase tol
kernel.alpha=1;%shape parameter. Small alpha-instability, large alpha-low accuracy
%divergence free matrix
D=points;%position of velocity vectors, a m-by-3 matrix
C=points;%position of Gaussian centers, a n-by-3 matrix
G=Gaussian_matrix(D,C,kernel);%for simplicity, D and C coincide
b=zeros(3*size(C,1),1);
b(1:3:end)=vel_in(:,1);
b(2:3:end)=vel_in(:,2);
b(3:3:end)=vel_in(:,3);

%% If C and D conincide, G is a symmetric and square matrix and can be solved by Conjugate Gradient method
% However if C and D do not coincide, G is a rectangular matrix and should
% be solved by GMRES or other methods
lamma=pcg(G,b,tol*stopping_constant);%conjugate gradient only works for symmetric matrix

%% Now place any points within the pipe to reconstruct its velocity, you can generate randomly
p2=[];
[x,y]=meshgrid(-cylinder.radius:x_w:cylinder.radius,-cylinder.radius:y_w:cylinder.radius);
idx=x.^2+y.^2<=cylinder.radius^2;
x=x(idx);
y=y(idx);
for z=cylinder.length_limit(1)+spacing+0.5:spacing:cylinder.length_limit(2)-spacing-0.5
    p2=[p2;x(:),y(:),z*ones(numel(x),1)];
end
v2=v_poiseulle(p2,cylinder,Q);
figure;
colored_quiver(p2(:,1),p2(:,2),p2(:,3),v2(:,1),v2(:,2),v2(:,3),scale,' ');
title('ground truth at the testing points');

%%
G2=Gaussian_matrix(p2,C,kernel);
v2_out=G2*lamma;
u=v2_out(1:3:end);
v=v2_out(2:3:end);
w=v2_out(3:3:end);
figure;
colored_quiver(p2(:,1),p2(:,2),p2(:,3),u,v,w,scale,' ');
title('reconstructed velocity at the testing points');

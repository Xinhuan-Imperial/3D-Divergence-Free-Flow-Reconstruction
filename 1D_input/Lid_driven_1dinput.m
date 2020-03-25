% by Dr. Xinhuan Zhou from Imperial College London, zhouxinhuan0205@126.com
% This code is for 3D reconstruction of steady Lid Driven Cavity flow
% using 1D velocity vector components. Cavity length:4mm-4mm-4mm. The
% lid(y=4) moves along z-axis at 0.15m/s.
% note:3D reconstruction from 1D velocity vectors seems only work when strong vortex exists. It needs to be
% tested for experimental feasibility. The algorithm is originally designed
% for cortex flow in ventricles. Generally speaking, as the input only has
% 1D component with minimal information, 3D reconstruction from 1D has
% limited accuracy - but this method can use a combination of 1D, 2D and 3D
% input.

clear;close all;clc;
scale=5;%scaling parameter for Matlab quiver function

%% extract uniform 3D velocity from non-uniform 3D lid driven cavity flow CFD simulation results.
% Lid driven cavity flow: top lid moving at 1m/s, CFD obtained by STAR-CCM+
load('LDC.mat');
%spacing in x,y and z-axis, uniform for simplicity but not necessary
x_w=0.4e-3;
y_w=0.4e-3;
z_w=0.4e-3;

%interpolate the CFD results to a uniform grid with grid spacing
%x_w,y_w,z_w
[X,Y,Z]=meshgrid(min(x):x_w:max(x),min(y):y_w:max(y),min(z):z_w:max(z));
X=X(:);Y=Y(:);Z=Z(:);
FU = scatteredInterpolant(x,y,z,u,'natural');%interpolate to uniform grid
FV = scatteredInterpolant(x,y,z,v,'natural');
FW = scatteredInterpolant(x,y,z,w,'natural');
U=FU(X,Y,Z);
V=FV(X,Y,Z);
W=FW(X,Y,Z);

h=figure(1);
set(h,'position',[0 0 1200 800]);
subplot(121);colored_quiver(x,y,z,u,v,w,scale,' ');daspect([1 1 1])
title('before downsampling');
subplot(122);colored_quiver(X,Y,Z,U,V,W,scale,' ');daspect([1 1 1])
title('after downsampling: ground truth');

%% assume the Doppler probe is at Probe. As the probe moves further, less information will be extracted so that the probe should be close to the flow field.
%find the points at the six walls
points2=[];
[x_tmp,y_tmp]=meshgrid(0:x_w:4e-3,0:y_w:4e-3);
points2=[points2;x_tmp(:),y_tmp(:),0*ones(numel(x_tmp),1)];
points2=[points2;x_tmp(:),y_tmp(:),4e-3*ones(numel(x_tmp),1)];
[x_tmp,z_tmp]=meshgrid(0:x_w:4e-3,0:z_w:4e-3);
points2=[points2;x_tmp(:),0*ones(numel(x_tmp),1),z_tmp(:)];
points2=[points2;x_tmp(:),4e-3*ones(numel(x_tmp),1),z_tmp(:)];
[y_tmp,z_tmp]=meshgrid(0:y_w:4e-3,0:z_w:4e-3);
points2=[points2;0*ones(numel(y_tmp),1),y_tmp(:),z_tmp(:)];
vel_in2=zeros(size(points2));
points2=[points2;4e-3*ones(numel(y_tmp),1),y_tmp(:),z_tmp(:)];%this wall has a velocity of 0.15m/s along a-axis
vel_in2=[vel_in2;zeros(numel(y_tmp),1),zeros(numel(y_tmp),1),0.15*ones(numel(y_tmp),1)];

%the position of probe is very important- must extract as much information
%as possible, and is flow dependent
Probe=[0,2e-3,2e-3];
points1=[X,Y,Z];%points1 are internal velocity vectors
vel_in1=[];
R=zeros(3*size(X,1)+3*size(points2,1),3*size(X,1)+3*size(points2,1));
for i=1:size(X,1)
    B=points1(i,:)-Probe;
    B=B/norm(B);%tmp is a unit vector pointing from Probe to the velocity vector
    tmp_vel=[U(i),V(i),W(i)];
    P=[B(1)*B(1),B(1)*B(2),B(1)*B(3);
        B(1)*B(2),B(2)*B(2),B(2)*B(3);
        B(1)*B(3),B(2)*B(3),B(3)*B(3)]; %projection matrix for tmp_vel
    tmp_vel=(P*tmp_vel')';%the projected velocity along the tmp_pos
    vel_in1=[vel_in1;tmp_vel];
    R(3*i-2:3*i,3*i-2:3*i)=P;
end

points=[points1;points2];
vel_in=[vel_in1;vel_in2];
m1=size(points1,1);
m2=size(points2,1);
R(3*m1+1:3*m1+3*m2,3*m1+1:3*m1+3*m2)=kron(eye(m2),eye(3,3));

figure;
plot3(Probe(1),Probe(2),Probe(3),'r*');hold on;%the red point is the artificial Doppler probe
colored_quiver(points(:,1),points(:,2),points(:,3),vel_in(:,1),vel_in(:,2),vel_in(:,3),scale,' ');daspect([1 1 1])
title('projected velocity:input 1D velocity vectors and 3D wall velocity');

%% construct the DFI matrix
stopping_constant=1.05;%early-stopping constant
%tol is very important to avoid overfitting, as input contains error
tol=5e-1;%input has large error due to: interpolation, CFD simulation
kernel.alpha=1/mean([x_w,y_w,z_w]);%shape parameter. Small alpha-instability, large alpha-low accuracy
%divergence free matrix
D=points;%position of velocity vectors, a m-by-3 matrix
C=points;%position of Gaussian centers, a n-by-3 matrix
G=Gaussian_matrix(D,C,kernel);%for simplicity, D and C coincide

%%
b=zeros(3*size(D,1),1);
b(1:3:end)=vel_in(:,1);
b(2:3:end)=vel_in(:,2);
b(3:3:end)=vel_in(:,3);
A=R*G;
max_it=m1+3*m2;
lamma=gmres(A,b,[],tol*stopping_constant,max_it);%conjugate gradient only works for symmetric matrix

%% Now place any points within the pipe to reconstruct its velocity, you can generate randomly
p2=[];
[x2,y2,z2]=meshgrid(min(x)+x_w/2:x_w:max(x)-x_w/2,min(y)+y_w/2:y_w:max(y)-y_w/2/2,min(z)+z_w/2:z_w:max(z)-z_w/2);
p2=[x2(:),y2(:),z2(:)];
G2=Gaussian_matrix(p2,C,kernel);
U2=FU(p2(:,1),p2(:,2),p2(:,3));
V2=FV(p2(:,1),p2(:,2),p2(:,3));
W2=FW(p2(:,1),p2(:,2),p2(:,3));

%%
v2_out=G2*lamma;
h=figure;
set(h,'position',[100 100 1000 600]);
subplot(121)
colored_quiver(p2(:,1),p2(:,2),p2(:,3),U2,V2,W2,scale,' ');
title('ground truth 3D velocity');
subplot(122)
colored_quiver(p2(:,1),p2(:,2),p2(:,3),v2_out(1:3:end),v2_out(2:3:end),v2_out(3:3:end),scale,' ');
title('reconstructed 3D velocity');


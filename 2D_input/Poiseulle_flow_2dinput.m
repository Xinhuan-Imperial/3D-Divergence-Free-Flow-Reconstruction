%% by Dr. Xinhuan Zhou from Imperial College London, zhouxinhuan0205@126.com
% This code is for 3D reconstruction of steady Poiseulle flow in a ciruclar
% pipe, using 2D velocity input
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
theta1=45/180*pi;%two imaging angles
theta2=-45/180*pi;

T1=[cos(theta1) 0 sin(theta1);0 1 0;-sin(theta1) 0 cos(theta1)];%rotation matrix
T2=[cos(theta2) 0 sin(theta2);0 1 0;-sin(theta2) 0 cos(theta2)];
B1=T1*[0;0;1];B2=-T2*[0;0;1];
skew_B1=[0,-B1(3),B1(2);B1(3),0,-B1(1);-B1(2),B1(1),0];
skew_B2=[0,-B2(3),B2(2);B2(3),0,-B2(1);-B2(2),B2(1),0];
R11=-skew_B1*skew_B1/single(norm(B1,2)^2);%projection matrix
R22=-skew_B2*skew_B2/single(norm(B2,2)^2);

%% Calculate input 2D velocity vectors: assume 2D velocities are planaer projection of 3D velocities
[x,y]=meshgrid(-cylinder.radius*3:x_w:cylinder.radius*3,-cylinder.radius*3:y_w:cylinder.radius*3);
points1=[x(:),y(:),zeros(numel(x),1)];
points1=(T1*points1')';
idx=points1(:,1).^2+points1(:,2).^2<=cylinder.radius^2;
points1=points1(idx,:);
points1(:,3)=points1(:,3)-max(points1(:,3))+cylinder.length_limit(1);

points2=[x(:),y(:),ones(numel(x),1)];
points2=(T2*points2')';
idx=points2(:,1).^2+points2(:,2).^2<=cylinder.radius^2;
points2=points2(idx,:);
points2(:,3)=points2(:,3)-max(points2(:,3))+cylinder.length_limit(1);

theta=linspace(0,2*pi,20)';
points3=[cylinder.radius*cos(theta),cylinder.radius*sin(theta),zeros(size(theta))];

points=[];
vel_in1=[];
vel_in2=[];
vel_in3=[];
% points from imaging angle1
for z=cylinder.length_limit(1):spacing:cylinder.length_limit(2)*2
    point_tmp1=points1;
    point_tmp1(:,3)=point_tmp1(:,3)+z;
    idx=point_tmp1(:,3)>=cylinder.length_limit(1) & point_tmp1(:,3)<=cylinder.length_limit(2);
    point_tmp1=point_tmp1(idx,:);
    points=[points;point_tmp1];  
    vel1=v_poiseulle(point_tmp1,cylinder,Q);
    vel1=(R11*vel1')';%conduct projection
    vel_in1=[vel_in1;vel1];
end

% points from imaging angle2
for z=cylinder.length_limit(1):spacing:cylinder.length_limit(2)*2
    point_tmp2=points2;
    point_tmp2(:,3)=point_tmp2(:,3)+z;
    idx=point_tmp2(:,3)>=cylinder.length_limit(1) & point_tmp2(:,3)<=cylinder.length_limit(2);
    point_tmp2=point_tmp2(idx,:);
    points=[points;point_tmp2];   
    vel2=v_poiseulle(point_tmp2,cylinder,Q);
    vel2=(R22*vel2')';    
    vel_in2=[vel_in2;vel2];
end

% points from surface scanning
for z=cylinder.length_limit(1):spacing:cylinder.length_limit(2)*2
    point_tmp3=points3;
    point_tmp3(:,3)=point_tmp3(:,3)+z;
    idx=point_tmp3(:,3)>=cylinder.length_limit(1) & point_tmp3(:,3)<=cylinder.length_limit(2);
    point_tmp3=point_tmp3(idx,:);
    points=[points;point_tmp3];     
    vel_in3=[vel_in3;zeros(size(point_tmp3))];
end

vel_in=[vel_in1;vel_in2;vel_in3];
m1=size(vel_in1,1);
m2=size(vel_in2,1);
m3=size(vel_in3,1);
figure;
colored_quiver(points(:,1),points(:,2),points(:,3),vel_in(:,1),vel_in(:,2),vel_in(:,3),scale,' ');
title('input 2D velocity- projected to two imaging planes');

%% construct the DFI matrix
stopping_constant=1.05;%early-stopping constant
tol=1e-2;%input is free of error, so tol can be very small here. If error exists, increase tol
kernel.alpha=1;%shape parameter. Small alpha-instability, large alpha-low accuracy
%divergence free matrix
D=points;%position of velocity vectors, a m-by-3 matrix
C=points;%position of Gaussian centers, a n-by-3 matrix
G=Gaussian_matrix(D,C,kernel);%for simplicity, D and C coincide
b=zeros(3*size(D,1),1);
b(1:3:end)=vel_in(:,1);
b(2:3:end)=vel_in(:,2);
b(3:3:end)=vel_in(:,3);
%NB: R is the projection matrix. Here its first 3*m1 diagonal elements corresponds to angle1, then 3m2 for angle2.
% the last 3*m3 elements correspond to surface points
R=zeros(3*m1+3*m2+3*m3,3*m1+3*m2+3*m3);
R(1:3*m1,1:3*m1)=kron(eye(m1),R11);
R(3*m1+1:3*m1+3*m2,3*m1+1:3*m1+3*m2)=kron(eye(m2),R22);
R(3*m1+3*m2+1:end,3*m1+3*m2+1:end)=kron(eye(m3),eye(3,3));

%% A=R*G is not symmetric, and should use GMRES with proper regularization
A=R*G;
lamma=gmres(A,b,[],tol*stopping_constant);

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

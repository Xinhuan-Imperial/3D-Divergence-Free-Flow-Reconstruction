% by Dr. Xinhuan Zhou from Imperial College London, zhouxinhuan0205@126.com
% This code is for 3D reconstruction of unsteady Womersley flow in a ciruclar pipe, using 3D velocity input
clear;close all;clc;

%% parameters of the straight pipe and input velocity vector distribution
cylinder.length_limit=[0 10];%total length along z-axis
cylinder.center=[0 0];
cylinder.radius=5;

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
vel_in=womersley_func(points(:,1)/1000,points(:,2)/1000,points(:,3)/1000);
h=figure;
set(h,'position',[100 100 1000 600]);
for i=1:size(vel_in,2)
    colored_quiver(points(:,1),points(:,2),points(:,3),zeros(size(points,1),1),zeros(size(points,1),1),vel_in(:,i),scale,' ');
    title(['unsteady input 3D velocity: ' num2str(i)]);
    pause(0.1);
end

%% construct the DFI matrix
stopping_constant=1.05;%early-stopping constant
tol=1e-3;%input is free of error, so tol can be very small here. If error exists, increase tol
kernel.alpha=1;%shape parameter. Small alpha-instability, large alpha-low accuracy
%divergence free matrix
D=points;%position of velocity vectors, a m-by-3 matrix
C=points;%position of Gaussian centers, a n-by-3 matrix
G=Gaussian_matrix(D,C,kernel);%for simplicity, D and C coincide

%% If C and D conincide, G is a symmetric and square matrix and can be solved by Conjugate Gradient method
% However if C and D do not coincide, G is a rectangular matrix and should
% be solved by GMRES or other methods
b=zeros(3*size(C,1),1);
for i=1:size(vel_in,2)
    b(1:3:end)=zeros(size(points,1),1);
    b(2:3:end)=zeros(size(points,1),1);
    b(3:3:end)=vel_in(:,i);
    lamma(:,i)=pcg(G,b,tol*stopping_constant);%conjugate gradient only works for symmetric matrix
end

%% Now place any points within the pipe to reconstruct its velocity, you can generate randomly
p2=[];
[x,y]=meshgrid(-cylinder.radius:x_w:cylinder.radius,-cylinder.radius:y_w:cylinder.radius);
idx=x.^2+y.^2<=cylinder.radius^2;
x=x(idx);
y=y(idx);
for z=cylinder.length_limit(1)+spacing+0.5:spacing:cylinder.length_limit(2)-spacing-0.5
    p2=[p2;x(:),y(:),z*ones(numel(x),1)];
end
v2=womersley_func(p2(:,1)/1000,p2(:,2)/1000,p2(:,3)/1000);

%%
G2=Gaussian_matrix(p2,C,kernel);
for i=1:size(vel_in,2)
    v2_out(:,i)=G2*lamma(:,i);
end

h=figure;
set(h,'position',[100 100 1000 600]);
for i=1:size(vel_in,2)
    subplot(121)
    colored_quiver(p2(:,1),p2(:,2),p2(:,3),zeros(size(p2,1),1),zeros(size(p2,1),1),v2(:,i),scale,' ');
    title(['unsteady ground truth 3D velocity: ' num2str(i)]);
    subplot(122)
    colored_quiver(p2(:,1),p2(:,2),p2(:,3),v2_out(1:3:end,i),v2_out(2:3:end,i),v2_out(3:3:end,i),scale,' ');
    title(['reconstructed 3D velocity: ' num2str(i)]);
    pause(0.1);
end

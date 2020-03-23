% by Dr. Xinhuan Zhou from Imperial College London, zhouxinhuan0205@126.com
% This code is for 3D reconstruction of unsteady Womersley flow in a
% ciruclar pipe, using 2D velocity vectors
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

points2=[x(:),y(:),zeros(numel(x),1)];
points2=(T2*points2')';
idx=points2(:,1).^2+points2(:,2).^2<=cylinder.radius^2;
points2=points2(idx,:);
points2(:,3)=points2(:,3)-max(points2(:,3))+cylinder.length_limit(1);

theta=linspace(0,2*pi,30)';
points3=[cylinder.radius*cos(theta),cylinder.radius*sin(theta),zeros(size(theta))];

points=[];
points_tmp=[];
% points from imaging angle1
for z=cylinder.length_limit(1):spacing:cylinder.length_limit(2)*2
    point_tmp1=points1;
    point_tmp1(:,3)=point_tmp1(:,3)+z;
    idx=point_tmp1(:,3)>=cylinder.length_limit(1) & point_tmp1(:,3)<=cylinder.length_limit(2);
    point_tmp1=point_tmp1(idx,:);
    points_tmp=[points_tmp;point_tmp1];
end
vel_in1=womersley_func(points_tmp(:,1)/1000,points_tmp(:,2)/1000,points_tmp(:,3)/1000);
points=[points;points_tmp];

points_tmp=[];
% points from imaging angle2
for z=cylinder.length_limit(1):spacing:cylinder.length_limit(2)*2
    point_tmp2=points2;
    point_tmp2(:,3)=point_tmp2(:,3)+z;
    idx=point_tmp2(:,3)>=cylinder.length_limit(1) & point_tmp2(:,3)<=cylinder.length_limit(2);
    point_tmp2=point_tmp2(idx,:);
    points_tmp=[points_tmp;point_tmp2];
end
vel_in2=womersley_func(points_tmp(:,1)/1000,points_tmp(:,2)/1000,points_tmp(:,3)/1000);
points=[points;points_tmp];

points_tmp=[];
% points from surface scanning, smaller spacing because velocity has large
% gradient near wall for Womersley flow
for z=cylinder.length_limit(1):spacing/5:cylinder.length_limit(2)*2
    point_tmp3=points3;
    point_tmp3(:,3)=point_tmp3(:,3)+z;
    idx=point_tmp3(:,3)>=cylinder.length_limit(1) & point_tmp3(:,3)<=cylinder.length_limit(2);
    point_tmp3=point_tmp3(idx,:);
    points_tmp=[points_tmp;point_tmp3];
end
points=[points;points_tmp];
vel_in3=zeros(size(points_tmp,1),size(vel_in1,2));

uin1=[];vin1=[];win1=[];uin2=[];vin2=[];win2=[];uin3=[];vin3=[];win3=[];
for i=1:size(vel_in1,2)
    tmp=[zeros(1,size(vel_in1,1));zeros(1,size(vel_in1,1));vel_in1(:,i)'];
    tmp=(R11*tmp)';
    uin1=[uin1,tmp(:,1)];
    vin1=[vin1,tmp(:,2)];
    win1=[win1,tmp(:,3)];
    tmp=[zeros(1,size(vel_in2,1));zeros(1,size(vel_in2,1));vel_in2(:,i)'];
    tmp=(R22*tmp)';
    uin2=[uin2,tmp(:,1)];
    vin2=[vin2,tmp(:,2)];
    win2=[win2,tmp(:,3)];   
    uin3=[uin3,zeros(size(vel_in3,1),1)];
    vin3=[vin3,zeros(size(vel_in3,1),1)];
    win3=[win3,zeros(size(vel_in3,1),1)];     
end
u_in=[uin1;uin2;uin3];
v_in=[vin1;vin2;vin3];
w_in=[win1;win2;win3];

m1=size(vel_in1,1);
m2=size(vel_in2,1);
m3=size(vel_in3,1);

h=figure;
set(h,'position',[100 100 1000 600]);
for i=1:size(u_in,2)
    colored_quiver(points(:,1),points(:,2),points(:,3),u_in(:,i),v_in(:,i),w_in(:,i),scale,' ');daspect([1 1 1])
    title(['unsteady input 2D velocity: ' num2str(i)]);
    pause(0.1);
end

%% construct the DFI matrix
stopping_constant=1.05;%early-stopping constant
tol=0.2;%input is the numerical solution of Womersley function and contains error due to truncation. To avoid overfitting, tol is 0.2
kernel.alpha=1;%shape parameter. Small alpha-instability, large alpha-low accuracy
%divergence free matrix
D=points;%position of velocity vectors, a m-by-3 matrix
C=points;%position of Gaussian centers, a n-by-3 matrix
G=Gaussian_matrix(D,C,kernel);%for simplicity, D and C coincide

R=zeros(3*m1+3*m2+3*m3,3*m1+3*m2+3*m3);
R(1:3*m1,1:3*m1)=kron(eye(m1),R11);
R(3*m1+1:3*m1+3*m2,3*m1+1:3*m1+3*m2)=kron(eye(m2),R22);
R(3*m1+3*m2+1:end,3*m1+3*m2+1:end)=kron(eye(m3),eye(3,3));

%% If C and D conincide, G is a symmetric and square matrix and can be solved by Conjugate Gradient method
% However if C and D do not coincide, G is a rectangular matrix and should
% be solved by GMRES or other methods
b=zeros(3*size(D,1),1);
lamma=zeros(3*size(C,1),size(u_in,2));
A=R*G;
max_it=2*m1+2*m2+3*m3;
for i=1:size(u_in,2)
    b(1:3:end)=u_in(:,i);
    b(2:3:end)=v_in(:,i);
    b(3:3:end)=w_in(:,i);
    lamma(:,i)=gmres(A,b,[],tol*stopping_constant,max_it);%conjugate gradient only works for symmetric matrix
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
G2=Gaussian_matrix(p2,C,kernel);

%%
for i=1:size(u_in,2)
    v2_out(:,i)=G2*lamma(:,i);
end

h=figure;
set(h,'position',[100 100 1000 600]);
for i=1:size(u_in,2)
    subplot(121)
    colored_quiver(p2(:,1),p2(:,2),p2(:,3),zeros(size(p2,1),1),zeros(size(p2,1),1),v2(:,i),scale,' ');
    title(['unsteady ground truth 3D velocity: ' num2str(i)]);
    subplot(122)
    colored_quiver(p2(:,1),p2(:,2),p2(:,3),v2_out(1:3:end,i),v2_out(2:3:end,i),v2_out(3:3:end,i),scale,' ');
    title(['reconstructed 3D velocity: ' num2str(i)]);
    pause(0.1);
end

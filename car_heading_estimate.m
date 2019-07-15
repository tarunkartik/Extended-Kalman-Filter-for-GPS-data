
close all
clear all
clc
load('gps_data.mat');
N1=200;
N=360;
t=data(N1+1:N1+N,1); % time 
%% true data from GPS
xtrue=data(N1+1:N1+N,2);
ytrue=data(N1+1:N1+N,3);
ztrue=data(N1+1:N1+N,4);
% sensor data with Noise
xsensor=xtrue+randn(1)*1;
ysensor=ytrue+randn(1)*1;
zsensor=ztrue+randn(1)*1;
%N=size(t,1);

pn = xtrue(1);
pe = ytrue(1);
xhat = [pn;pe;15;15];
measure = [xsensor ysensor];
%% states xfilter=[x,y,vx,vy]
% Filter Initialization
% xhat = [xsensor; ysensor; 0;0];
P = [...
    5^2 0 0 0;...
    0 5^2 0 0 ;...
    0 0 7^2 0 ;...
    0 0 0 7^2];
Q = 1*eye(4);
rn = 0.2;
re = 0.2;
R = [rn 0;
    0 re];

A = [...
    0 0 1 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];

B = [...
    0 0;
    0 0;
    1 0;
    0 1];

C = [1 0 0 0;
    0 1 0 0];
u = [0,0;
    0, 0];

var =xhat'; 
     for i=2:N
         dt=t(i)-t(i-1);% compute time diff between two nodes
         M = 10;
         %dt=0.1;
         for j = 1:M
         %% prediction
        xhat = xhat + (dt/M)*(A*xhat+B*u);
        P = P + (dt/M)*(A*P + P*A' + Q);
        P_track(:,:,i)= P;
        covariance_pn(i) = P_track(1,1,i);
        covariance_pe(i) = P_track(2,2,i);
        covariance_vn(i) = P_track(3,3,i);
%         covariance_ve(i) = P_track(4,4,i);
        
         end
         %% Update
%          cn = C(1,:);
%          ce = C(2,:);
         L = P*C'*(R+C*P*C')^(-1);
%          Le = P*ce'*(re+ce*P*ce')^(-1);
%          xhat = xhat+ Ln*();
         P = (eye(4)-L*C)*P;
         xhat = xhat + L*(measure(i,:)'-xhat(1:2,1));
         
         var = vertcat(var,xhat(:,1)');
     end
     
  error_pn = xtrue(1:360,1)- var(1:360,1);
  error_pe = ytrue(1:360,1)- var(1:360,2);
%   error_vn = vntrue(1:N,1)- var(:,1);
%   error_pn = pntrue(1:N,1)- var(:,1);

ubound_pn = 3*sqrt(covariance_pn);
lbound_pn = -3*sqrt(covariance_pn);
  

ubound_pe = 3*sqrt(covariance_pe);
lbound_pe = -3*sqrt(covariance_pe);
  
  figure(1)
  plot(xtrue, ytrue,'-b')
  hold on
  grid on
  plot(var(:,1),var(:,2), '--r')
  title('Trajectory Plot')
%   xlabel('time (in sec)')
%   ylabel('error in pn')
  legend('True Position', 'Estimated Position')
  
  figure(2)
  plot(ubound_pn, '-r')
  hold on
  grid on
  plot(error_pn, 'b')
  hold on
  plot(lbound_pn, 'r' )
  title('3-sigma Plot for pn')
  xlabel('time (in sec)')
  ylabel('error in pn')
  legend('Bounds', 'Error')
  
  figure(3)
  plot(ubound_pe, '-r')
  hold on
  grid on
  plot(error_pe, 'b')
  hold on
  plot(lbound_pe, 'r' )
  title('3-sigma Plot for pe')
  xlabel('time (in sec)')
  ylabel('error in pe')
  legend('Bounds', 'Error')
  
  figure(4)
  plot(xtrue(:,1),'-r')
  hold on
  grid on
  plot(var(:,1),'-b')
  title('True v/s estimated position for pn')
  xlabel('time (in sec)')
  ylabel('Car position in north')
  legend('True position', 'Estimated position')
  
  figure(5)
  plot(ytrue(:,1),'-r')
  hold on
  grid on
  plot(var(:,2),'-b')
  title('True v/s estimated position for pe')
  xlabel('time (in sec)')
  ylabel('Car position in east')
  legend('True position', 'Estimated position')
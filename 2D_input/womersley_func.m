% generate the Womersley flow based on some preset parameters by Xiaowei
% Zhou 2018.July

function Results=womersley_func(x,y,z)
rho=1.06e3;          %  Density of blood   [kg/m^3]
mu=0.0035;            %  Viscosity of blood [kg/m s]
R=3/1000;          %  Radius of vessel   [m]
Heartbeat = 50;   % the heartbeat in beats/min
omega_0=2*pi*Heartbeat/60;  %  Angular frequency
V0=1/15;             %  Mean velocity      [m/s]. Note: this mean velocity is the temporal mean velocity during one cardiac cycle which is calculated from the pulsatile spatial mean velocity curve
r=sqrt(x.^2+y.^2);        %  Radial sampling interval 0.002, when use this script for generating the reference velocity, deltaR = 0.02;
PRF = Heartbeat;
r_rel=r/R;  %  Values for the relative radius
numsample = length(r_rel);
%   Calculate psi
for p=1:8
    
  disp(['Harmonic number ',num2str(p)])
  omega=p*omega_0;
  tau_alpha=1i^(3/2)*R*sqrt(omega*rho/mu);  
  Be=tau_alpha*bes(0,tau_alpha);  % bes is the same as the MATLAB besselj function, but bes was customly written.
  psi(:,p)=(Be-tau_alpha*bes(0,r_rel*tau_alpha))/(Be-2*bes(1,tau_alpha));
  
end

% the following Vp is the mean velocity waveform information in the Fourier
% domain. The column 1 is the amplitude, the column 2 is the phase. The
% first 8 components of the waveform in the Fourier Domain were used.
Vp=[1.00    0.00         % the zeroth frequency magnitude is 1, which means the spatial pulsatile mean velocity is 1 m/s 
    1.89   32
    2.49   85
    1.28  156
    0.32  193
    0.27  133
    0.32  155
    0.28  195
    0.01  310];

Vp(:,1) = Vp(:,1)*V0;  % in this way, the mean velocity V0 whose value is kind of a percentage value (or scale factor) for the 1 m/s in Vp
Vp(:,2) = Vp(:,2)*pi/180;
time=linspace(0,60/60,PRF);  % only one cycle of the data will be simulated, which is 1s.
len=length(time);

%  Generate files for the scatteres over a number of pulse emissions
Nshoots=len;
F=struct('cdata',[],'colormap',[]);
for i=1:Nshoots
    %  Calculate the profile
    prof=2*Vp(1,1)*(1-r_rel.^2); % prof is the velocity along the radius from centre to vessel wall
    index=1:numsample;
    for p=2:9
        prof=prof + Vp(p,1)*abs(psi(index,p-1)).*cos((p-1)*omega_0*time(i) - Vp(p,2) + angle(psi(index,p-1)));
    end
%     plot(prof);drawnow;
    Results(:,i) = prof;   % Results is the profile for the whole cardiac cycle
end




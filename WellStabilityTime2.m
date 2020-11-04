%% Azad_Almasov

%% Advanced_Drilling_-_2018_-_Fall

% Final Project

clc
clear all

%%#######################################################################%%
% TIME DELAYED WELLBORE FAILURE:
%% Given Data:

D         =  9482;      % ft, true vertical depth
phi_deg   =  42.4;      % degree of angle of internal friction
A_deg     =  10;        % degree of angle of azimuth
grad      =  1;         % psi/ft, vertical stress gradient

biot     =  1;         % Biott Coefficient
S_0      =  1182;      % psi, cohesive strength
Pi       =  3200;      % psi, pore pressure
nu_h     =  0.25;     % horizontal Poisson's ratio
nu_vmin  =  0.2;      % minimum vertical Poisson's ratio
nu_vmax  =  0.35;      % maximum vertical Poisson's ratio

T_0      =  0;         % psi, tensile strength
P_w      =  5500;      % psi, wellbore pressure

inc_deg   = (0.1:5:90.1);  % degree of inlination angle
theta_deg = 90; % degree of radial change

%%
Cs = 1.81e-7;   % Solid compessibility, 1/psi
Cb = 8.62e-7;   % Bulk compressibility, 1/psi
Cw = 2.86e-6;   % Water compressibility, 1/psi
qw = - 100000;    % Water flow rate, STB/D
h = 160;        % Resservoir thickness, ft
rr = 327;       % Observation well distance, ft
phi = 0.145;    % Porosity
Bw = 1;         % Water formation volume factor Bbl/STB
mu = 1.16;      % Viscosity of water, cP
k =  75;        % Permeability, mD;
Case      =  input('Cases: isotropic(1) and anisotropic(2)    ');
case1     =  input('Cases for BCs: plain strain (1); generalized plain stress     ');
%%#######################################################################%%
%% SOLUTION

t = logspace(log10(0.01),log10(50),5)*365*24; % Production time in hrs.
% Plain strain storage capacity:

if case1 ==1      % (plane strain)
S = (-2/3)*(((Cb-Cs)^2)/Cb)*((1-2*nu_h)/(1-nu_h))+phi*Cw-(1+phi)*Cs+Cb;
else if case1==2  % (generalized plane stress)
S = ((Cb-Cs)/3)*((1-Cs/Cb)*((3-4*nu_h)*(1+nu_h)/(2*(1-nu_h^2))-3)) ...
    +Cb-Cs*(1+phi)+phi*Cw;
    end
end
zeta = 0.00264*k/(mu*S);

% Degrees to Radians
phi_rad   = phi_deg*pi/180;
A_rad     = A_deg*pi/180;
inc_rad   = inc_deg*pi/180;
theta_rad = theta_deg*pi/180;
% Calculation of sigma_v:
sigma_v   = grad*D;

for tt=1:length(t);
% Pore pressure changes with time. Line source solution:
P_pore(tt) = Pi - (70.6*qw*mu/(k*h))*expint((rr^2)/(4*zeta*t(tt)));

% Stress path:
if Case   == 1          % Isotropic case    
sigma_h0  =  6143;      % psi
sigma_h(tt)   =  sigma_h0 - biot*(1-2*nu_h)*(P_pore(tt) -Pi)/(1-nu_h);
sigma_H(tt)   =  sigma_h(tt);

else                    % Anisotropic case              
sigma_h0   =  6143;
sigma_h(tt)   =  sigma_h0 - biot*(1-(nu_vmin+nu_h*nu_vmax)/(1-nu_h^2))*(P_pore(tt) -Pi);

sigma_H0   =  7500;
sigma_H(tt)   =  sigma_H0 - biot*(1-(nu_vmin*nu_h+nu_vmax)/(1-nu_h^2))*(P_pore(tt) -Pi);
end

%% Transformation to Wellbore cartesian system:
for j=1:length(inc_deg)

    
sigma_x{tt}(j)  = (sigma_H(tt)*(cos(A_rad))^2 + sigma_h(tt)*(sin(A_rad))^2) ... 
              *(cos(inc_rad(j)))^2 + sigma_v*(sin(inc_rad(j)))^2;

sigma_y{tt}(j)  = sigma_H(tt)*(sin(A_rad))^2 + sigma_h(tt)*(cos(A_rad))^2;

sigma_z{tt}(j)  = (sigma_H(tt)*(cos(A_rad))^2 + sigma_h(tt)*(sin(A_rad))^2) ... 
              *(sin(inc_rad(j)))^2 + sigma_v*(cos(inc_rad(j)))^2;

sigma_xy{tt}(j) = 0.5*(sigma_h(tt) - sigma_H(tt))*sin(2*A_rad)*cos(inc_rad(j));

sigma_xz{tt}(j) = 0.5*(sigma_H(tt)*(cos(A_rad))^2 + sigma_h(tt)*(sin(A_rad))^2 ...
              - sigma_v)*sin(2*inc_rad(j));

sigma_yz{tt}(j) = 0.5*(sigma_h(tt) - sigma_H(tt))*sin(2*A_rad)*sin(inc_rad(j));


%#########################################################################%
%% Fracture:
theta_rad_min{tt}(j) = 0.5*atan(2*sigma_xy{tt}(j)/(sigma_x{tt}(j) - sigma_y{tt}(j)));

sigma_tzf{tt}(j)     = 2*(-sigma_xz{tt}(j)*sin(theta_rad_min{tt}(j)) + sigma_yz{tt}(j) ...
                   *cos(theta_rad_min{tt}(j)));

P_wf{tt}(j)          = sigma_x{tt}(j) + sigma_y{tt}(j) - 2*(sigma_x{tt}(j) ...
                   - sigma_y{tt}(j))*cos(2*theta_rad_min{tt}(j)) - ...
                   4*sigma_xy{tt}(j)*sin(2*theta_rad_min{tt}(j)) ...
                   - (sigma_tzf{tt}(j)^2)/(sigma_z{tt}(j) - T_0 - biot*P_pore(tt)) ...
                   -biot*P_pore(tt) - T_0;

ppg_max{tt}(j)       = P_wf{tt}(j)/(0.052*D);



%#########################################################################%
%% Transform to radial coordinate system:
  

sigma_rr{tt}(j)  = P_w;    

sigma_zz{tt}(j) = sigma_z{tt}(j) -2*nu_h*(sigma_x{tt}(j) - sigma_y{tt}(j))*cos(2 ...
                 *theta_rad) - 4*nu_h*sigma_xy{tt}(j)*sin(2*theta_rad);

sigma_tt{tt}(j) = (sigma_x{tt}(j) + sigma_y{tt}(j)) - 2*(sigma_x{tt}(j) - sigma_y{tt}(j))...
                  *cos(2*theta_rad) - 4*sigma_xy{tt}(j)*sin(2* ...
                  theta_rad)-P_w;

sigma_tz{tt}(j) = 2*(-sigma_xz{tt}(j)*sin(theta_rad) + sigma_yz{tt}(j) ...
                 *cos(theta_rad));

sigma_rt{tt}(j) = 0;
 
sigma_rz{tt}(j) = 0;



%#########################################################################%
%% Collapse:
sigma_1c{tt}(j) = 0.5*(sigma_tt{tt}(j) + sigma_zz{tt}(j)) + 0.5 ...
                 *sqrt((sigma_tt{tt}(j) - sigma_zz{tt}(j))^2 ...
                 + 4*(sigma_tz{tt}(j))^2) - biot*P_pore(tt);

sigma_3c{tt}(j) = P_w - biot*P_pore(tt);

% Mohr Circle:
sigma_1{tt}(j) = 0.5*(sigma_tt{tt}(j) + sigma_zz{tt}(j)) + 0.5 ...
                 *sqrt((sigma_tt{tt}(j) - sigma_zz{tt}(j))^2 ...
                 + 4*(sigma_tz{tt}(j))^2) - biot*P_pore(tt);
             
sigma_3{tt}(j) = 0.5*(sigma_tt{tt}(j) + sigma_zz{tt}(j)) - 0.5 ...
                 *sqrt((sigma_tt{tt}(j) - sigma_zz{tt}(j))^2 ...
                 + 4*(sigma_tz{tt}(j))^2) - biot*P_pore(tt);
             
R{tt}(j) = (sigma_1c{tt}(j)-sigma_3c{tt}(j))/2;
end

center_ver{tt} = (sigma_1c{tt}(1)+sigma_3c{tt}(1))/2;
viscircles([center_ver{tt},0],R{tt}(1));
title('Mohr circle in vertical well');
xlabel('\sigma , psi');
ylabel('\tau , psi');
hold on

% center_inc{tt} = (sigma_1c{tt}(4)+sigma_3c{tt}(4))/2;
% viscircles([center_inc{tt},0],R{tt}(4));
% title('Mohr circle in vertical well');
% xlabel('\sigma , psi');
% ylabel('\tau , psi');
% hold on

% Coulomb plot
sigma_11(tt) = sigma_1c{tt}(19);
sigma_33(tt) = sigma_3c{tt}(19);
end

figure
plot(sigma_33,sigma_11);
xlabel('\sigma_3, psi');
ylabel('\sigma_1, psi');

% % Mohr circle for different inclination angles:
%  for tt=1:length(t);
% figure
%      for i=1:2:7
% center_inc{tt} = 0.5*(sigma_tt{tt}(i)+sigma_zz{tt}(i))-P_pore(tt);
% viscircles([center_inc{tt},0],R{tt}(i));
% title('Mohr circle for different inclination angles wells');
% xlabel('\sigma , psi');
% ylabel('\tau , psi');
% hold on
%      end
%  end
   
 ty = t/365/24;                   % Time in years.
 
 for i=1:4;
     for tt=1:length(t);
 Rt{i}(tt) = R{tt}(i);
     end
figure
 plot(ty, Rt{i})
 xlabel('time, years')
 ylabel('\tau_{max}, psi')
 hold on
 end

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%#############################___THE_END___#########################%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
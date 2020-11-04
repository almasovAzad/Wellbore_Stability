%% Azad_Almasov

%% Advanced_Drilling_-_2018_-_Fall

% Final Project

clc
clear all

%%#######################################################################%%
% STATEMENT OF THE PROBLEM 
%% Given Data:

D         =  9482;      % ft, true vertical depth
phi_deg   =  42.4;      % degree of angle of internal friction
A_deg     =  10;        % degree of angle of azimuth
grad      =  1;         % psi/ft, vertical stress gradient

biot     =  1;         % Biott Coefficient
S_0      =  1182;      % psi, cohesive strength
Pi       =  3500;      % psi, pore pressure
nu_h     =  0.25;     % horizontal Poisson's ratio
nu_vmin  =  0.2;      % minimum vertical Poisson's ratio
nu_vmax  =  0.35;      % maximum vertical Poisson's ratio

T_0      =  0;         % psi, tensile strength
syms     P_w;          % psi, wellbore pressure

inc_deg   = (0.1:5:90.1);  % degree of inlination angle
theta_deg = 90; % degree of radial change

%%
Cs = 1.81e-7;   % Solid compessibility, 1/psi
Cb = 8.62e-7;   % Bulk compressibility, 1/psi
Cw = 2.86e-6;   % Water compressibility, 1/psi
qw = 2734;      % Water flow rate, STB/D
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



% Mohr-Coulomb equation:
f   = inline(S_0+(0.5*(sigma_1c{tt}(j)+sigma_3c{tt}(j))-0.5 ...
             *(sigma_1c{tt}(j) - sigma_3c{tt}(j))*sin(phi_rad)) ...
             *tan(phi_rad) - 0.5*(sigma_1c{tt}(j) - sigma_3c{tt}(j)) ...
             *cos(phi_rad),'P_w');


%#########################################################################%
%% Newton_-_Raphson_-_Iteration:
df   = inline(diff(f(P_w)),'P_w');
ddf  = inline(diff(df(P_w)),'P_w');

acc  = 0.001;
P_wi = Pi-1500;

y1   = f(P_wi);
y2   = df(P_wi);
y3   = ddf(P_wi);
   
% Constion for the convertion criteria:  
   while (((y1*y3)/(y2*y2))>1)
      P_wi = input('Enter initial value for Pressure again')
      y1   = f(P_wi);
      y2   = df(P_wi);
      y3   = ddf(P_wi);
   end
P_w2{tt}(j) = (P_wi - (y1/y2));
% To check for the accuracy (convergence)
    while (abs(P_w2{tt}(j) - P_wi)>acc)
      P_wi = P_w2{tt}(j);
      y1   = f(P_wi);
      y2   = df(P_wi);
     
P_w2{tt}(j) =(P_wi - (y1/y2));
    end
    


P_wc{tt}(j) = P_w2{tt}(j);
ppg_min{tt}(j) = P_wc{tt}(j)/(0.052*D);

end


%#########################################################################%
%% Plot


plot(inc_deg,ppg_min{tt});
hold on

plot(inc_deg,ppg_max{tt});
hold on

if Case == 1
    title('Mud Window for the isotropic case (Case 1)');
else
    title('Mud Window for the anisotropic case (Case 2)');
end
xlabel('Inclination Degree');
ylabel('Mud Density');
end


ty = t/365/24;       % Time in years
figure
plot(ty, P_pore)
xlabel('time, years');
ylabel('P_{pore} , psi');

figure
plot(ty, sigma_H)
xlabel('time, years');
ylabel('\sigma_H , psi');

figure
plot(ty, sigma_h)
xlabel('time, years');
ylabel('\sigma_h , psi');

for kk=1:length(t)
ppg_max_vertical(kk) = ppg_max{kk}(1);
ppg_min_vertical(kk) = ppg_min{kk}(1);
dppg_vertical(kk) = ppg_max{kk}(1)- ppg_min{kk}(1);

ppg_max_horizontal(kk) = ppg_max{kk}(length(inc_deg)-1);
ppg_min_horizontal(kk) = ppg_min{kk}(length(inc_deg)-1);
dppg_horizontal(kk) = ppg_max{kk}(length(inc_deg)-1)- ppg_min{kk}(length(inc_deg)-1);
end

figure
plot(ty, dppg_vertical);
xlabel('time, years');
ylabel('mud weight difference in vertical well, ppg');

figure
plot(ty, dppg_horizontal);
xlabel('time, years')
ylabel('mud weight difference in horizontal well, ppg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%#############################___THE_END___#########################%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
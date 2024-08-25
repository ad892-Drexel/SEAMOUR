%% clear

close all
clc
tic
%% Trajectory Sender

%run("Trajectory_Sender_20221028.m")
addpath('Simscape_Parts_StepFiles\')
%% Important Initial Inputs

period=3;
simulation_ts=0.001;
speed_limit=600;
Tf=9;

%% MAIN BODY
test=0.2;
testz = 0;
gravity = 9.81;
testvel=0;

hinge_stiffness=.00001;
hinge_damp=.03;

% hinge_stiffness=10;
% hinge_damp=0.1;
offset=0.002; % bouyancy

%%
Ix_hexcore = 0.000392;
Iy_hexcore = 1.1548103e-01;
Iz_hexcore = 1.1548103e-01;
mass_hexcore = 3.48;

Ix_flooded =  6.8508989e-02;
Iy_flooded = 4.2964749e-01 ;
Iz_flooded = 4.2964749e-01 ;
mass_flooded = 1.3488164e+01;

dist_x = 0.062;
dist_y = 0;
dist_z = testz;

newIx_hexcore = Ix_hexcore+(mass_hexcore*(dist_x^2));
newIy_hexcore = Iy_hexcore+(mass_hexcore*(dist_y^2));
newIz_hexcore = Iz_hexcore+(mass_hexcore*(dist_z^2));

newIx_flooded = Ix_flooded+(mass_flooded*(dist_x^2));
newIy_flooded = Iy_flooded+(mass_flooded*(dist_y^2));
newIz_flooded = Iz_flooded+(mass_flooded*(dist_z^2));

final_Ix = newIx_hexcore + newIx_flooded;
final_Iy = newIy_hexcore + newIy_flooded;
final_Iz = newIz_hexcore + newIz_flooded;

center_body.mass = mass_flooded+mass_hexcore;
center_body.com = [0 0 0];
center_body.inertia = [final_Ix final_Iy final_Iz];
center_body.prodinertia = [0 0 0];

Body_Length=50;

%% Head

headmount.mass = 0;
headmount.com = [0 0 0];
headmount.inertia = [0 0 0];
headmount.prodinertia = [0 0 0];

headmotorshield.mass = 0;
headmotorshield.com = [0 0 0];
headmotorshield.inertia = [0 0 0];
headmotorshield.prodinertia = [0 0 0];

head_Ix_flooded = 4.0866377e-03;
head_Iy_flooded = 4.7789219e-03;
head_Iz_flooded = 4.7789219e-03;
head_mass_flooded = 1.7934437e+00 ;

head_Ix = (0.715*3.076e-3)/0.5;
head_Iy = (0.715*3.365e-3)/0.5;
head_Iz = (0.715*3.365e-3)/0.5;
head_mass = 0.715;

total_Ix_head = head_Ix_flooded+head_Ix;
total_Iy_head = head_Iy_flooded+head_Iy;
total_Iz_head = head_Iz_flooded+head_Iz;

headmain.mass = head_mass+head_mass_flooded;
headmain.com = [0 0 0];
headmain.inertia = [total_Ix_head, total_Iy_head, total_Iz_head];
headmain.prodinertia = [0 0 0];

%% Hind

tailmount.mass = 0;
tailmount.com = [0 0 0];
tailmount.inertia = [0 0 0];
tailmount.prodinertia = [0 0 0];

tailmotorshield.mass = 0;
tailmotorshield.com = [0 0 0];
tailmotorshield.inertia = [0 0 0];
tailmotorshield.prodinertia = [0 0 0];


Ix_tail_flooded = 7.1327049e-03;
Iy_tail_flooded = 1.1448029e-02;
Iz_tail_flooded = 1.1448029e-02;
mass_tail_flooded =  2.7530755e+00;

tail_Ix = (1.31*3.076e-3);
tail_Iy = (1.31*3.365e-3);
tail_Iz = (1.31*3.365e-3);
mass_tail = 1.31;

total_Ix_tail = Ix_tail_flooded + tail_Ix;
total_Iy_tail = Iy_tail_flooded + tail_Iy;
total_Iz_tail = Iz_tail_flooded + tail_Iz;


hindconnection.mass = mass_tail+mass_tail_flooded;
hindconnection.com = [-.110 0 0];
hindconnection.inertia = [total_Ix_tail, total_Iy_tail, total_Iz_tail];
hindconnection.prodinertia = [0 0 0];

hindcover.mass = 0;
hindcover.com = [0 0 0];
hindcover.inertia = [0 0 0];
hindcover.prodinertia = [0 0 0];


%% Left Flipper (27cmx7cmx2cm)

leftmotor1.mass = 0;
leftmotor1.com = [0 0 0];
leftmotor1.inertia = [0 0 0];
leftmotor1.prodinertia = [0 0 0];

leftmotor2.mass = 0;
leftmotor2.com = [0 0 0];
leftmotor2.inertia = [0 0 0];
leftmotor2.prodinertia = [0 0 0];


leftmotor3.mass = 0;
leftmotor3.com = [0 0 0];
leftmotor3.inertia = [0 0 0];
leftmotor3.prodinertia = [0 0 0];


lf_lower.mass = 0.072;
lf_lower.com = [0 0 0];
lf_lower.moi=[6.15e-5 2.20e-5 7.95e-5];
lf_lower.poi=[3.4e-7 0.111e-7 4.69e-7];

lf_hinge.stiffness = hinge_stiffness;
lf_hinge.damping = hinge_damp;

lf_upper.mass = 0.096;
lf_upper.com = [0 0 0];
lf_upper.moi=[18.1e-5 2.88e-5 20.8e-5];
lf_upper.poi=[-2.68e-10 1.58e-11 9.904e-6];


%% Right Flipper

rightmotor3.mass = 0;
rightmotor3.com = [0 0 0];
rightmotor3.inertia = [0 0 0];
rightmotor3.prodinertia = [0 0 0];

rightmotor2.mass = 0;
rightmotor2.com = [0 0 0];
rightmotor2.inertia = [0 0 0];
rightmotor2.prodinertia = [0 0 0];

rightmotor1.mass = 0;
rightmotor1.com = [0 0 0];
rightmotor1.inertia = [0 0 0];
rightmotor1.prodinertia = [0 0 0];

rf_lower.mass = 0.072;
rf_lower.com = [0 0 0];
rf_lower.moi=[6.15e-5 2.20e-5 7.95e-5];
rf_lower.poi=[-3.4e-7 0.111e-7 -4.69e-7];

rf_hinge.stiffness = hinge_stiffness;
rf_hinge.damping = hinge_damp;

rf_upper.mass = 0.096;
rf_upper.com = [0 0 0];
rf_upper.moi=[18.1e-5 2.88e-5 20.8e-5];
rf_upper.poi=[-2.68e-10 1.58e-11 9.904e-6];


%% Left Hind Flipper(18cmx7cmx1cm)

lefthindmotor.mass = 0;
lefthindmotor.com = [0 0 0];
lefthindmotor.inertia = [0 0 0];
lefthindmotor.prodinertia = [0 0 0];

lhfcomx = 0.112179;
lhfcomy = 0.0133764;
lhfcomz = -0.0128043;
lhf_inertiax = 6.03879e-05;
lhf_inertiay = 0.000404492;
lhf_inertiaz = 0.000348104;
lhf_prodinertiax = -1.65586e-08;
lhf_prodinertiay =-5.02123e-06;
lhf_prodinertiaz = -1.19244e-06;

lefthindflipper.mass = 0.132;
lefthindflipper.com = [0 0 0];
lefthindflipper.inertia = [lhf_inertiax, lhf_inertiay ,lhf_inertiaz];
lefthindflipper.prodinertia = [lhf_prodinertiax ,lhf_prodinertiay, lhf_prodinertiaz];

%% Right Hind Flipper

righthindmotor.mass = 0;
righthindmotor.com = [0 0 0];
righthindmotor.inertia = [0 0 0];
righthindmotor.prodinertia = [0 0 0];

rhfcomx = 0.112179;
rhfcomy = 0.0133764;
rhfcomz = -0.0128043;
rhf_inertiax = 6.03879e-05;
rhf_inertiay = 0.000404492;
rhf_inertiaz = 0.000348104;
rhf_prodinertiax = -1.65586e-08;
rhf_prodinertiay =-5.02123e-06;
rhf_prodinertiaz = -1.19244e-06;

righthindflipper.mass = 0.132;
righthindflipper.com = [0 0 0];
righthindflipper.inertia = [rhf_inertiax rhf_inertiay rhf_inertiaz];
righthindflipper.prodinertia = [rhf_prodinertiax rhf_prodinertiay rhf_prodinertiaz];




%% ADDED MASS FOR ELLIPSOID
a=0.25;
b=0.125;
rho=997;
m = (4*pi*rho*a*b*b)/3;
ecc = sqrt(1-((b/a)^2));
alpha0 = ((2*(1-(ecc^2)))/(ecc^3))*((0.5*log((1+ecc)/(1-ecc))-ecc));
beta0 = (1/(ecc^2))-(((1-(ecc^2))/(2*(ecc^3)))*(log((1+ecc)/(1-ecc))));
k1 = (alpha0)/(2-alpha0);
k2 = (beta0)/(2-beta0);
xudot = k1*m;
yvdot = k2*m;
zwdot = k2*m;

%% ADDED MASS FOR FOREFLIPPERS

length_lower= 0.09;
length_upper= 0.18;
thickness= 0.07;
width= 0.022;
[kx_upper,ky_upper,kz_upper]=AM_Flat_Plate(length_upper,thickness,width);
[kx_lower,ky_lower,kz_lower]=AM_Flat_Plate(length_lower,thickness,width);

%% ******************* Calculate Addeded Mass Distributions***************
% Main Body Surface Area
x_rad=.490;
y_rad=.125;
z_rad=.125;
sections=10;
nums=1000;

%% Send added mass distrbutions to simscape

% This array the first point is head and the last point is hind
Cd_y = [0.5645,0.5232,0.4085,0.3589,0.3589,0.40849,0.5232,0.4642];
Cd_z = [0.5645,0.52,0.5843,0.533,0.533,0.5843,0.52,0.5238];
Area_y = [0.0031581,0.024875,0.025026,0.021123,0.021123,0.025026,0.024875,0.037486];
Area_z = [0.0031588,0.028932,0.035488,0.030129,0.030129,0.035488,0.028932,0.039984];
t_y = sum(Area_y);
percent_y=Area_y/t_y;
t_z = sum(Area_z);
percent_z=Area_z/t_z;
Coeff=[xudot, yvdot, zwdot];



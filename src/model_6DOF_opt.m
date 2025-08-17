clc
close all
clear variables


g = 9.81;

lp = 0.85; % vzdalenost T - predni naprava
lz = 0.7; % vzdalenost T - zadni naprava
rp = 0.5; % predni rozchod/2
rz = 0.6; % zadni rozchod/2
h = 130e-3; % vyska ulozeni pruzeni v systemu sasi

R = 0.2341; % polomer kola
f0 = 0.7;

iw = 0.2372; % moment setrvacnosti kol
Ix = 20.22; % moment setrvacnosti chassis k x
Iy = 68.81; % moment setrvacnosti chassis k y
Iz = 71.42; % moment setrvacnosti chassis k z
Dxz = 5.8713; % deviacni moment xz

m = 320;
mw = 12.5;

M = diag([m; m; m]);
I = [ Ix, 0, -Dxz; 0, Iy, 0; -Dxz, 0, Iz];
Iw = diag([iw; iw; iw; iw]);

kp = 20000;
kz = 20000;
bp = 200;
bz = 200;

z0 = 300e-3; % vyska teziste od zeme

%% integrace

F_prl_p = m*g*lz/(2*(lp+lz));
F_prl_z = m*g/2 - F_prl_p;

finalt = 75;

dt = 1e-4;
t = 0 : 1e-4 : finalt;

ode_set = odeset('RelTol',1e-3, 'AbsTol',1e-4, 'maxstep', 1e-3);
    
options = optimset(	'LargeScale', 'off', ...
                   	'Display', 'iter', ...
                   	'TolFun', 1e-6, ...
                  	'TolX', 1e-6,...
                   	'MaxFunEvals', 30000,...
                    'MaxIter', 30000);
delta_max = 10;
M_max = 21;

R_des = 10;
v_des = 30;

t_init = 50;

% -------------------
% FWD, RWD
op0 = [ 8.8/delta_max; 
        10/M_max;
        15/M_max];
    
op_LB = zeros(3,1);                                                  
op_UB = ones(3,1);

% -------------------

% pro fixni natoceni kol
% op_LB(1) = 1; 

% -------------------

u0 = [  v_des/3.6;0;0;
        0;0;0;
        1e-3;1e-3;1e-3;1e-3;...
        0;0;300e-3;
        0;0;0;
        0;0;0;0]; %poc. podminky
    
save('data.mat')  
                
[op_py, fval, exitflag, output, lambda] = fmincon(@objfce_opt, op0, [], [], [], [], op_LB, op_UB, [], options);   
        
op_py = op_py.*[delta_max; M_max; M_max];

function objfcn = objfce_opt(op)
                
load('data.mat')

delta = delta_max*op(1);

% RWD
M_1 = 0;
M_2 = 0;
M_3 = M_max*op(2);
M_4 = M_max*op(3);

% % FWD
% M_1 = M_max*op(2);
% M_2 = M_max*op(3);
% M_3 = 0;
% M_4 = 0;


[T, U] = ode15s(@(t,u)model_6DoF(t, u, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0, delta, M_1, M_2, M_3, M_4), t, u0, ode_set);

ind_init = t_init/dt;

stred_X = mean(U(ind_init : end, 11));
stred_Y = mean(U(ind_init : end, 12));

R_akt = sqrt((U(ind_init : end,11) - stred_X).^2 + (U(ind_init : end,12) - stred_Y).^2);

objfcn = max((1 - R_akt/R_des).^2);

% objfcn = max((1 - U(ind_init : end,4)/1.17).^2);

% objfcn =    max((1 - R_akt/R_des).^2) + ...
%             max((1 - U(ind_init : end,4)/(v_des/3.6/R_des)).^2);

end


function du = model_6DoF(t, u, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0, delta, M_1, M_2, M_3, M_4)

dxdt = u(1);
dydt = u(2);
dzdt = u(3);
dpsidt = u(4);
dthetadt = u(5);
dfidt = u(6);
omega1 = u(7);
omega2 = u(8);
omega3 = u(9);
omega4 = u(10);
x = u(11);
y = u(12);
z = u(13);
psi = u(14);
theta = u(15);
fi = u(16);
fi1 = u(17);
fi2 = u(18);
fi3 = u(19);
fi4 = u(20);

%% transformacni matice z lokalniho systemu do globalniho

% yaw angle
A_psi = [   cos(psi), -sin(psi), 0; 
            sin(psi), cos(psi), 0;
            0, 0, 1];
     
% pitch angle
A_theta = [ cos(theta), 0, sin(theta); 
            0, 1, 0;
            -sin(theta), 0, cos(theta)];
       
% roll angle
A_fi = [    1, 0, 0; 
            0, cos(fi), -sin(fi); 
            0, sin(fi), cos(fi)];
        
A = A_psi*A_theta*A_fi;
A_theta_fi = A_theta*A_fi;

% matice G        
G = [   0, -sin(psi), cos(theta)*cos(psi);
        0, cos(psi), cos(theta)*sin(psi);
        1, 0, -sin(theta)];
    
dGdt = [    0, -dpsidt*cos(psi), (-dthetadt*sin(theta))*cos(psi)+cos(theta)*(-dpsidt*sin(psi));
            0, -dpsidt*sin(psi), (-dthetadt*sin(theta))*sin(psi)+cos(theta)*(dpsidt*cos(psi));
            0, 0, -dthetadt*cos(theta)];
        
% polohove vektory uchyceni pruzeni   
r_1 = [lp; rp; h];
r_2 = [lp; -rp; h];
r_3 = [-lz; rz; h];
r_4 = [-lz; -rz; h];

% deformace odpruzeni 
def_1 = A_theta_fi*r_1 - r_1;
def_2 = A_theta_fi*r_2 - r_2;
def_3 = A_theta_fi*r_3 - r_3;
def_4 = A_theta_fi*r_4 - r_4;

% rychlost deformace odpruzeni 
def_1dt = u(1:3) + cross(G*u(4:6), A*r_1);
def_2dt = u(1:3) + cross(G*u(4:6), A*r_2);
def_3dt = u(1:3) + cross(G*u(4:6), A*r_3);
def_4dt = u(1:3) + cross(G*u(4:6), A*r_4);

r_1k = [r_1(1:2) + def_1(1:2); -z];
r_2k = [r_2(1:2) + def_2(1:2); -z];
r_3k = [r_3(1:2) + def_3(1:2); -z];
r_4k = [r_4(1:2) + def_4(1:2); -z];

% sily v odpruzeni
F_s1 = F_prl_p - kp*(z - z0 + def_1(3)) - bp*def_1dt(3); 
F_s2 = F_prl_p - kp*(z - z0 + def_2(3)) - bp*def_2dt(3); 
F_s3 = F_prl_z - kz*(z - z0 + def_3(3)) - bz*def_3dt(3); 
F_s4 = F_prl_z - kz*(z - z0 + def_4(3)) - bz*def_4dt(3);

F_s = [F_s1; F_s2; F_s3; F_s4];

N_k = mw*g + F_s;

%% rizeni momentu

% [M_1, M_2, M_3, M_4] = rizeni_momentu(t);

%% natoceni kol   

delta1 = deg2rad(delta);
delta2 = deg2rad(delta);

%% speed, slip angles, tire force

ro_f = atan(rp/lp);
ro_r = atan(rz/lz);

v = sqrt(dxdt^2 + dydt^2);
d_f = sqrt(rp^2 + lp^2);
d_r = sqrt(rz^2 + lz^2);

Beta = atan2(dydt,dxdt) - psi;

v_lo(1) = (-dpsidt*d_f*sin(ro_f) + v*cos(Beta));
v_lo(2) = ( dpsidt*d_f*sin(ro_f) + v*cos(Beta));
v_lo(3) = (-dpsidt*d_r*sin(ro_r) + v*cos(Beta));
v_lo(4) = ( dpsidt*d_r*sin(ro_r) + v*cos(Beta));

v_tr(1) = ( dpsidt*d_f*cos(ro_f) + v*sin(Beta));
v_tr(2) = ( dpsidt*d_f*cos(ro_f) + v*sin(Beta));
v_tr(3) = (-dpsidt*d_r*cos(ro_r) + v*sin(Beta));
v_tr(4) = (-dpsidt*d_r*cos(ro_r) + v*sin(Beta));

v_wh_ax(1) = v_lo(1)*cos(delta1) + v_tr(1)*sin(delta1);
v_wh_ax(2) = v_lo(2)*cos(delta2) + v_tr(2)*sin(delta2);
v_wh_ax(3) = v_lo(3);
v_wh_ax(4) = v_lo(4);

beta(1) = delta1 - atan2(v_tr(1),v_lo(1));
beta(2) = delta2 - atan2(v_tr(2),v_lo(2));
beta(3) = - atan2(v_tr(3),v_lo(3));
beta(4) = - atan2(v_tr(4),v_lo(4));

C_sigma = 50000;
C_S = 35000;
zz = 0.025;

F_z = N_k;

for i = 1 : 4
    
    denom = max(abs(R*u(6+i)),abs(v_wh_ax(i)));
    sigma(i) = (R*u(6+i) - v_wh_ax(i))/denom;

    lambda(i) = f0*F_z(i)*(1+sigma(i))/(2*sqrt((C_sigma*sigma(i))^2+(C_S*tan(beta(i)))^2));
    if (lambda(i) < 1)
         f_lambda(i) = (2-lambda(i))*lambda(i);
    else
         f_lambda(i) = 1;
    end

    H(i) = C_sigma*f_lambda(i)*sigma(i)/(1+sigma(i));
    S(i) = C_S*f_lambda(i)*tan(beta(i))/(1+sigma(i));
    Mk(i) = zz*S(i);

end

%% rhs

% sily pusobici na ram
F_1 = [cos(delta1), -sin(delta1), 0; sin(delta1), cos(delta1), 0; 0, 0, 1]*[H(1); S(1); F_s(1)];
F_2 = [cos(delta2), -sin(delta2), 0; sin(delta2), cos(delta2), 0; 0, 0, 1]*[H(2); S(2); F_s(2)];
F_3 = [H(3); S(3); F_s(3)];
F_4 = [H(4); S(4); F_s(4)];

f = A_psi*(F_1 + F_2 + F_3 + F_4) - [0; 0; m*g];

% G_arb = 8.1e10;
% a_arb = 100e-3;
% l_arb = 500e-3;
% D_arb = 16e-3;
% d_arb = 14e-3;
% J_arb = (pi*(D_arb^4 - d_arb^4))/32;
% k_arb = G_arb*J_arb*l_arb/(a_arb^2);

% momenty sil pusobici na ram
fm =    cross(A_theta_fi'*r_1k, A_theta_fi'*F_1) + ...
        cross(A_theta_fi'*r_2k, A_theta_fi'*F_2) + ...
        cross(A_theta_fi'*r_3k, A_theta_fi'*F_3) + ...
        cross(A_theta_fi'*r_4k, A_theta_fi'*F_4) - ...
        A_theta_fi'*[0; (M_1 + M_2 + M_3 + M_4); 0] - ...
        A_theta_fi'*[0; 0; sum(Mk)];% ...
%         - [2*k_arb*fi; 0; 0];

omegapruh = A'*G*u(4:6);
omegaIomega = cross(omegapruh,(I*omegapruh));

% momenty sil na kolech
Mw = [  M_1 - R*H(1); 
        M_2 - R*H(2); 
        M_3 - R*H(3); 
        M_4 - R*H(4)];

% res(1,1) = u(1) - up(11);
% res(2,1) = u(2) - up(12);
% res(3,1) = u(3) - up(13);
% res(4,1) = u(4) - up(14);
% res(5,1) = u(5) - up(15);
% res(6,1) = u(6) - up(16);
% res(7,1) = u(7) - up(17);
% res(8,1) = u(8) - up(18);
% res(9,1) = u(9) - up(19);
% res(10,1) = u(10) - up(20);
% 
% res(11:13,1) = f - M*up(1:3);
% res(14:16,1) = fm - omegaIomega - I*(A'*(dGdt*u(4:6) + G*up(4:6)));
% res(17:20,1) = Mw - Iw*up(7:10);

du(11:20,1) = u(1:10,1);

du(1:3,1) = M\f;
du(4:6,1) = (I*A'*G)\(fm - omegaIomega - I*A'*dGdt*u(4:6,1));
du(7:10,1) = Iw\Mw;

end
    

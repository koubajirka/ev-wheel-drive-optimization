clc
clear variables
close all
%exportgraphics(figure(13), 'rizeni_rchlosti.pdf');
g = 9.81;

lp = 0.85; % vzdalenost T - predni naprava
lz = 0.7; % vzdalenost T - zadni naprava
rp = 0.5; % predni rozchod/2
rz = 0.6; % zadni rozchod/2
h = 130e-3; % vyska ulozeni pruzeni v systemu sasi

R = 0.2341; % polomer kola
f0 = 0.8;

vdes = 25;

iw = 0.2372; % moment setrvacnosti kol
Ix = 20.22; % moment setrvacnosti chassis k x
Iy = 68.81; % moment setrvacnosti chassis k y
Iz = 71.42; % moment setrvacnosti chassis k z
Dxz = 5.8713; % deviacni moment xz

m = 320;
mw = 12.5;

M = diag([m+4*mw; m+4*mw; m]);
I = [ Ix, 0, -Dxz; 0, Iy, 0; -Dxz, 0, Iz];
Iw = diag([iw; iw; iw; iw]);

kp = 10000;
kz = 10000;
bp = 1500;
bz = 1500;

z0 = 300e-3; % vyska teziste od zeme
%% modelo = 1 => rovne, modelo = 2 => kruznice
modelo = 2;


%% integrace

F_prl_p = m*g*lz/(2*(lp+lz));
F_prl_z = m*g/2 - F_prl_p;
finalt = 100;
if modelo == 1
finalt = 9;
end
dt = 1e-4;
t = 0 : 1e-4 : finalt;

ode_set = odeset('RelTol',1e-3, 'AbsTol',1e-4, 'maxstep', 1e-3);

u0 = [  vdes/3.6;0;0;
        0;0;0;
        1e-3;1e-3;1e-3;1e-3;...
        0;0;300e-3;
        0;0;0;
        0;0;0;0]; %poc. podminky
    
% up0 = zeros(20,1);
% fixed_u0 = ones(20,1);
% fixed_up0 = zeros(20,1);
% [u0_new, up0_new] = decic(@(t,u,up)model_6DoF(t, u, up, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0), t, u0, fixed_u0, up0, fixed_up0, ode_set);
% [T, U] = ode15i(@(t,u,up)model_6DoF(t, u, up, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0), t, u0, up0_new, ode_set);

[T, U] = ode15s(@(t,u)model_6DoF(t, u, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0,modelo), t, u0, ode_set);


%% vypocet relativnich rychlosti, slip angle etc

for ii = 1 : length(T)
    [delta1(ii,1), delta2(ii,1)] = natoceni_kol(T(ii),modelo);
end

ro_f = atan(rp/lp);
ro_r = atan(rz/lz);

v = sqrt(U(:,1).^2 + U(:,2).^2);
d_f = sqrt(rp^2 + lp^2);
d_r = sqrt(rz^2 + lz^2);

beta_auto = atan2(U(:,2),U(:,1)) - U(:,14);

v_lo(:,1) = (-U(:,4).*d_f*sin(ro_f) + v.*cos(beta_auto));
v_lo(:,2) = ( U(:,4).*d_f*sin(ro_f) + v.*cos(beta_auto));
v_lo(:,3) = (-U(:,4).*d_r*sin(ro_r) + v.*cos(beta_auto));
v_lo(:,4) = ( U(:,4).*d_r*sin(ro_r) + v.*cos(beta_auto));

v_tr(:,1) = ( U(:,4).*d_f*cos(ro_f) + v.*sin(beta_auto));
v_tr(:,2) = ( U(:,4).*d_f*cos(ro_f) + v.*sin(beta_auto));
v_tr(:,3) = (-U(:,4).*d_r*cos(ro_r) + v.*sin(beta_auto));
v_tr(:,4) = (-U(:,4).*d_r*cos(ro_r) + v.*sin(beta_auto));

v_wh_ax(:,1) = v_lo(:,1).*cos(delta1(:,1)) + v_tr(:,1).*sin(delta1(:,1));
v_wh_ax(:,2) = v_lo(:,2).*cos(delta2(:,1)) + v_tr(:,2).*sin(delta2(:,1));
v_wh_ax(:,3) = v_lo(:,3);
v_wh_ax(:,4) = v_lo(:,4);

beta(:,1) = delta1 - atan2(v_tr(:,1),v_lo(:,1));
beta(:,2) = delta2 - atan2(v_tr(:,2),v_lo(:,2));
beta(:,3) = - atan2(v_tr(:,3),v_lo(:,3));
beta(:,4) = - atan2(v_tr(:,4),v_lo(:,4));

skluz(:,1) = (R.*U(:,7) - v_wh_ax(:,1))./(max([R.*U(:,7),v_wh_ax(:,1)],[],2));
skluz(:,2) = (R.*U(:,8) - v_wh_ax(:,2))./(max([R.*U(:,8),v_wh_ax(:,2)],[],2));
skluz(:,3) = (R.*U(:,9) - v_wh_ax(:,3))./(max([R.*U(:,9),v_wh_ax(:,3)],[],2));
skluz(:,4) = (R.*U(:,10) - v_wh_ax(:,4))./(max([R.*U(:,10),v_wh_ax(:,4)],[],2));

if modelo == 1
 subplot(4,2,1)
set(gca, 'Fontsize', 10, 'layer', 'top', 'TickLength', [0.02,0.025])
hold on 
box on
plot(T, U(:,11), 'linewidth', 2)
xlabel('cas {\itt} (s)')
ylabel('vzdalenost {\itx} (m)')
title('draha')


subplot(4,2,2)
set(gca, 'Fontsize', 10, 'layer', 'top', 'TickLength', [0.02,0.025])
hold on 
box on
plot(T, U(:,1), 'linewidth', 2)
xlabel('cas {\itt} (s)')
ylabel('rychlost {\itv} (m/s)')
title('rychlost')



subplot(4,1,2)
set(gca, 'Fontsize', 10, 'layer', 'top', 'TickLength', [0.02,0.025])
hold on 
box on
plot(T, U(:,7), 'linewidth', 2)
plot(T, U(:,8), 'linewidth', 2, 'linestyle', '--')
plot(T, U(:,9), 'linewidth', 2, 'linestyle', '-.');
plot(T, U(:,10), 'linewidth', 2, 'linestyle', ':')
legend('omega PP','omega LP','omega PZ','omega LZ', 'Location', 'northeast')
xlabel('cas {\itt} (s)')
ylabel('uhlova rychlost {\it\omega} (rad/s)')
title('uhlove rychlosti')



subplot(4,2,5)
set(gca,  'Fontsize', 10, 'layer', 'top', 'TickLength', [0.02,0.025])
hold on 
box on
plot(T, U(:,5), 'linewidth', 2)
xlabel('cas {\itt} (s)')
ylabel('rychlost kloneni {\it\theta} (rad/s)')
title('rychlost kloneni {\it\theta}')


subplot(4,2,6)
set(gca,  'Fontsize', 10, 'layer', 'top', 'TickLength', [0.02,0.025])
hold on 
box on
plot(T, U(:,15), 'linewidth', 2)
xlabel('cas {\itt} (s)')
ylabel('uhel {\it\theta} (rad)')
title('kloneni')


subplot(4,1,4)
set(gca,  'Fontsize', 10, 'layer', 'top', 'TickLength', [0.02,0.025])
hold on 
box on
plot(T,R.*U(:,7)-U(:,1),'linewidth', 2);
plot(T,R.*U(:,8)-U(:,1),'linewidth', 2,'linestyle', '--');
plot(T,R.*U(:,9)-U(:,1),'linewidth', 2, 'linestyle', '-.');
plot(T,R.*U(:,10)-U(:,1),'linewidth', 2, 'linestyle', ':');
legend('kolo PP', 'kolo LP', 'kolo PZ', 'kolo LZ', 'Location', 'northeast')
xlabel('cas [s]')
ylabel('relativni rychlost [m/s]')   
end    
%%
if modelo == 2   
figure
hold on
plot(T,U(:,11), 'linewidth', 1);
plot(T,U(:,12), 'linewidth', 1);
plot(T,U(:,13), 'linewidth', 1);
ylabel('poloha vozu X,Y,Z [m]')
xlabel('cas [s]')
legend('X', 'Y', 'Z')

figure
hold on
plot(U(:,11),U(:,12),'linewidth', 1);
ylabel('poloha vozu Y [m]')
xlabel('poloha vozu X [m]')
axis equal

figure
hold on
plot(T,rad2deg(U(:,14)),'linewidth', 1);
legend('yaw')
ylabel('yaw [deg]')

figure
hold on
plot(T,rad2deg(U(:,15)),'linewidth', 1);
legend('pitch')
ylabel('pitch [deg]')

figure
hold on
plot(T,rad2deg(U(:,16)),'linewidth', 1);
legend('roll')
ylabel('roll [deg]')

figure
hold on
plot(T,3.6*U(:,1),'linewidth', 1);
plot(T,3.6*U(:,2),'linewidth', 1);
plot(T,3.6*U(:,3),'linewidth', 1);
ylabel('rychlosti vozu XY [km/h]')
xlabel('cas [s]')
legend('vX', 'vY', 'vZ')

figure
hold on
plot(T,U(:,4),'linewidth', 1);
plot(T,U(:,5),'linewidth', 1);
plot(T,U(:,6),'linewidth', 1);
legend('yaw rate', 'pitch rate' ,'roll rate')
ylabel('roll rate, pitch rate a yaw rate [rad/s]')

figure
hold on
plot(T,U(:,7),'linewidth', 1);
plot(T,U(:,8),'linewidth', 1);
plot(T,U(:,9),'linewidth', 1);
plot(T,U(:,10),'linewidth', 1);
legend('kolo 1', 'kolo 2', 'kolo 3', 'kolo 4')
xlabel('cas [s]')
ylabel('rychlost otaceni kol [rad/s]')

figure
hold on
plot(T,v_wh_ax(:,1),'linewidth', 1);
plot(T,v_wh_ax(:,2),'linewidth', 1);
plot(T,v_wh_ax(:,3),'linewidth', 1);
plot(T,v_wh_ax(:,4),'linewidth', 1);
legend('kolo 1', 'kolo 2', 'kolo 3', 'kolo 4')
xlabel('cas [s]')
ylabel('v_{wh_ax} [m/s]')

figure
hold on
plot(T,rad2deg(beta_auto(:,1)),'linewidth', 1);
xlabel('cas [s]')
ylabel('slip angle auto [deg]')

figure
hold on
plot(T,rad2deg(beta(:,1)),'linewidth', 1);
plot(T,rad2deg(beta(:,2)),'linewidth', 1);
plot(T,rad2deg(beta(:,3)),'linewidth', 1);
plot(T,rad2deg(beta(:,4)),'linewidth', 1);
legend('kolo 1', 'kolo 2', 'kolo 3', 'kolo 4')
xlabel('cas [s]')
ylabel('slip angle [deg]')

figure
hold on
plot(T,skluz(:,1),'linewidth', 1);
plot(T,skluz(:,2),'linewidth', 1);
plot(T,skluz(:,3),'linewidth', 1);
plot(T,skluz(:,4),'linewidth', 1);
legend('kolo 1', 'kolo 2', 'kolo 3', 'kolo 4')
xlabel('cas [s]')
ylabel('slip ratio [-]')

figure
plot(T, 3.6*sqrt(U(:,1).^2+U(:,2).^2))
xlabel('cas [s]')
ylabel('rychlost vozidla [km/h]')

ind_init = 35/dt;
stred_X = mean(U(ind_init : end, 11));
stred_Y = mean(U(ind_init : end, 12));
R_akt = sqrt((U(ind_init : end,11) - stred_X).^2 + (U(ind_init : end,12) - stred_Y).^2);
R_des = 10;

figure
plot(t(ind_init : end), R_akt-R_des)
xlabel('cas [s]')
ylabel('\Delta R [m]')

end

% function res = model_6DoF(t, u, up, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0)
function du = model_6DoF(t, u, M, I, Iw, m, mw, lp, lz, rp, rz, h, R, kp, kz, bp, bz, z0, F_prl_p, F_prl_z, g, f0,modelo)

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

[M_1, M_2, M_3, M_4] = rizeni_momentu(t,modelo);

%% natoceni kol   

[delta1, delta2] = natoceni_kol(t,modelo);

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

F_z(1:4) = N_k;

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
    

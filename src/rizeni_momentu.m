
function [M_1, M_2, M_3, M_4] = rizeni_momentu(t,modelo)
if modelo == 2
M_1 = 0; % moment na levem prednim kole
M_2 = 0; % moment na pravem prednim kole
M_3 = 8.6398; % moment na levem zadnim kole
M_4 = 12.6864; % moment na pravem zadnim kole
%--------------------1
M_3 = 5.9969;
M_4 = 6.1767;
%--------------------2
% M_3 = 5.8304;
% M_4 = 6.3267; %větší
%--------------------3
% M_3 = 8.2363;
% M_4 = 16.6097; %větší

% M_1 = 10; % moment na levem prednim kole
% M_2 = 10; % moment na pravem prednim kole
% M_3 = 10;
% M_4 = 10;
end
%------------------------------------------
if modelo == 1
M_1 = 45;
M_2 = 45;
M_3 = 45;
M_4 = 45;


if t>2
    M_1=M_1-M_1*(0.4*(t-2));
    M_2=M_2-M_2*(0.4*(t-2));
    M_3=M_3-M_3*(0.4*(t-2));
    M_4=M_4-M_4*(0.4*(t-2));
end

if  t>7
    M_1=-45;
    M_2=-45;
    M_3=-45;
    M_4=-45;
end
end
%------------------------------------------------
% zmeneny moment
% dM_1z = 0; % moment na levem prednim kole
% dM_2z = 0; % moment na pravem prednim kole
% dM_3z = 0; % moment na levem zadnim kole
% dM_4z = 0; % moment na pravem zadnim kole

% cas_zmeny = 30;
% delka_zmeny = 1; % casova delka zmeny momentu
% 
% if t > cas_zmeny && t <= (cas_zmeny + delka_zmeny)
% 
%     M_1 = M_1 + dM_1z*((t-cas_zmeny)/delka_zmeny).^2.*(3-2*((t-cas_zmeny)/delka_zmeny));   
%     M_2 = M_2 + dM_2z*((t-cas_zmeny)/delka_zmeny).^2.*(3-2*((t-cas_zmeny)/delka_zmeny));  
%     M_3 = M_3 + dM_3z*((t-cas_zmeny)/delka_zmeny).^2.*(3-2*((t-cas_zmeny)/delka_zmeny));  
%     M_4 = M_4 + dM_4z*((t-cas_zmeny)/delka_zmeny).^2.*(3-2*((t-cas_zmeny)/delka_zmeny));  
%     
% elseif t > (cas_zmeny + delka_zmeny)
%     
%     M_1 = M_1 + dM_1z;   
%     M_2 = M_2 + dM_2z;  
%     M_3 = M_3 + dM_3z;  
%     M_4 = M_4 + dM_4z;  
%     
% end

end
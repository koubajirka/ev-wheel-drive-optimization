
function [delta1, delta2] = natoceni_kol(t,modelo)
if modelo == 2;
uhel_1 = 8.5945;
uhel_1 = 8.7202;
%uhel_1 = 8.7164;
% uhel_1 = 8.9901;
uhel_2 = uhel_1;
end
if modelo == 1 
uhel_1 = 0;
uhel_2 = uhel_1;
end
delta1 = deg2rad(uhel_1);
delta2 = deg2rad(uhel_2);  

% duhel = 0;
% 
% cas_zatoceni = 20;
% delka_manevru = 0.5; % casova delka manevru

% if t > cas_zatoceni && t <= (cas_zatoceni+delka_manevru)
%     
%     delta1 = deg2rad(uhel_1) + deg2rad(duhel*((t-cas_zatoceni)/delka_manevru).^2.*(3-2*((t-cas_zatoceni)/delka_manevru)));   
%     delta2 = deg2rad(uhel_2);

% elseif t > (cas_zatoceni+delka_manevru)
%     
%     delta1 = deg2rad(uhel_1 + duhel);
%     delta2 = deg2rad(delta1);
%     
% end
    
end
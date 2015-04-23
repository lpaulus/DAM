function devoir

gamma = 1.3;
D = 0.075; % diametre piston [m]
R = 0.03675; % longueur manivelle [m]
L = 0.1285; % longueur bielle [m]
V_c = R*pi*D^2/2 % cylindree [m^3] (= 324.25 cm^3)
beta = L/R;
tau = 8.2; % taux de compression
Q_tot = 2800000; % chaleur degagee par la combustion [J/kg]
dTheta_comb = 40*pi/180; % duree combustion [rad]
theta_d = -15*pi/180; % debut combustion [rad]
w = [3000*2*pi/60 5000*2*pi/60]; % vitesses angulaire moteur [rad/s]
m_piston = 0.310; % [kg]
m_bielle = 0.544; % [kg]
p_atm = 101325; % [Pa]
p_0 = 10^5; % [bar]
rho_air = p_atm/(298.15*(8.3145/0.02896)); % [kg/m^3]
x_air = 14.7/15.7; % proportion d'air dans le cylindre
k = 1.5; % pour l'echappement
alpha = pi/2; % debut du temre de relaxation

    function dpdtheta = pression(theta, p) % fonction qui sera intégrée par ode45
        V = 0.5*V_c*(1-cos(theta)+beta-sqrt(beta^2 - (sin(theta))^2)+ (2/(tau-1)));
        dVdtheta = 0.5*V_c*(sin(theta) + (sin(theta)*cos(theta))/sqrt(beta^2 - (sin(theta))^2));
        
        if theta > theta_d && theta < theta_d + dTheta_comb % combustion
            dQdtheta = pi*Q_tot*V_c*tau/(tau-1)*x_air*rho_air*sin(pi*((theta-theta_d)/dTheta_comb))/(2*dTheta_comb);
        else
            dQdtheta = 0;
        end
        
        if theta <  -pi  % admission
            dpdtheta = 0;
        elseif theta > alpha  % echappement (équation non intégrée)
            dpdtheta = -k*(p-p_atm);
        else % compression et explosion
            dpdtheta = (gamma-1)*dQdtheta/V - gamma*p/V*dVdtheta;
        end
    end

options = odeset('MaxStep',1e-2);
[Theta, P] = ode45(@pression,[-2*pi;2*pi], p_atm, options);

%pour l'échappement version intégrée, même résultat
%{
for i=1:length(Theta)
   if Theta(i)>pi
       p_loc = P(i);
       i_loc = i;
       break;
   end
end

for j = i_loc:length(Theta) % echappement
   P(j) = p_atm+(p_loc-p_atm)*exp(-k*(Theta(j)-pi));
end
%}

Pp = P/p_0; % pour remettre en bars
p_max = max(Pp)
for i = 1:2 % forces
    F_pied(:,i) = pi*D^2*P/4 - m_piston*R*w(i)^2*cos(Theta);
    F_tete(:,i) = -pi*D^2*P/4 + (m_piston+m_bielle)*R*w(i)^2*cos(Theta);
    F_tot(:,i) = F_pied(:,i) - F_tete(:,i);
    
    F_max_min = [max(F_tete(:,i)) min(F_tete(:,i)) max(F_pied(:,i)) min(F_pied(:,i))]
    F_tot_max_min = [max(F_tot(:,i)) min(F_tot(:,i))]
    
    %{
    for j = 1:length(F_pied(:,i))
        if F_pied(j,i) < 0
            F_pied_trac(j,i) = - F_pied(j,i);
            F_pied_comp(j,i) = 0;
        else
            F_pied_trac(j,i) = 0;
            F_pied_comp(j,i) = F_pied(j,i);
        end
        if F_tete(j,i) < 0
            F_tete_comp(j,i) = - F_tete(j,i);
            F_tete_trac(j,i) = 0;
        else
            F_tete_comp(j,i) = 0;
            F_tete_trac(j,i) = F_tete(j,i);
        end
    end
    F_tot_comp(:,i) = F_pied_comp(:,i) + F_tete_comp(:,i);
    F_tot_trac(:,i) = F_pied_trac(:,i) + F_tete_trac(:,i);
    %}
end

% graphiques
close all;


% pression
figure;
plot(Theta, Pp);
title('Pression en fonction de l''angle du vilebrequin ')
xlabel('\theta [rad]')
ylabel('p [bar]')
% axes x et y
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','LineStyle',':','linewidth',0.2) %x-axis
line([0 0], yL,'color','k','LineStyle',':','linewidth',0.2) %y-axis
saveas(gcf, 'pression', 'epsc');

% efforts 3000 rpm
figure;
plot(Theta,F_pied(:,1)/1000,'-r',Theta,F_tete(:,1)/1000,'-b');
title({'Efforts sur la bielle en fonction de l''angle du vilebrequin';'f = 3000 rpm'})
legend('Pied','Tête')
xlabel('\theta [rad]')
ylabel('F [kN]')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','LineStyle',':','linewidth',0.2) %x-axis
line([0 0], yL,'color','k','LineStyle',':','linewidth',0.2) %y-axis
saveas(gcf, 'forces_3000rpm', 'epsc');

% effort total 3000 rpm
figure;
plot(Theta,F_tot(:,1)/1000);
title({'Efforts totaux sur la bielle en fonction de l''angle du vilebrequin';'f = 3000 rpm'})
xlabel('\theta [rad]')
ylabel('F [kN]')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','LineStyle',':','linewidth',0.2) %x-axis
line([0 0], yL,'color','k','LineStyle',':','linewidth',0.2) %y-axis
saveas(gcf, 'forces_tot_3000rpm', 'epsc');

% efforts 5000 rpm
figure;
plot(Theta,F_pied(:,2)/1000,'-r',Theta,F_tete(:,2)/1000,'-b');
title({'Efforts sur la bielle en fonction de l''angle du vilebrequin';'f = 5000 rpm'})
legend('Pied','Tête')
xlabel('\theta [rad]')
ylabel('F [kN]')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','LineStyle',':','linewidth',0.2) %x-axis
line([0 0], yL,'color','k','LineStyle',':','linewidth',0.2) %y-axis
saveas(gcf, 'forces_5000rpm', 'epsc');

% effort total 5000 rpm
figure;
plot(Theta,F_tot(:,2)/1000);
title({'Efforts totaux sur la bielle en fonction de l''angle du vilebrequin';'f = 5000 rpm'})
xlabel('\theta [rad]')
ylabel('F [kN]')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','LineStyle',':','linewidth',0.2) %x-axis
line([0 0], yL,'color','k','LineStyle',':','linewidth',0.2) %y-axis
saveas(gcf, 'forces_tot_5000rpm', 'epsc');

end
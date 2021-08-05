function [C_mean_major, C_mean_REE, P_mean, F_mean, crust_mean] = computeAverageProperties(M_out,F_out,C_fluid,C_fluid_major)
P = M_out(:,1);
F = F_out(:,1);
% H = grid_z(1)-grid_z(end);
% U = 2*((grid_z-grid_z(end))/H) - ((grid_z-grid_z(end))/H).^2;
  

% % Numerical solution for averaged melt composition - e.g. Ito 2005
% % Eq. 1
% % weighted with F
% FC = C_total_fluid_major.*F.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
% int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)));
% int_F = sum(((F(2:end).*U(2:end)+F(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)));
% C_ave_major_F = int_FC/int_F;

% mean pressure - Asimow 2001
int_P = cumsum(P.*[0 ;abs(F(2:end)-F(1:end-1))]); 
P_extr = 1./F.*int_P; P_extr(1) = 0;

FP = F.*P_extr; FP(1,:) = 0; FP(isnan(FP)) = 0;
int_FP = sum(((FP(2:end,:)+FP(1:end-1,:))/2).*abs(P(2:end)-P(1:end-1)));
int_F  = sum(((F(2:end)+F(1:end-1))/2).*abs(P(2:end)-P(1:end-1)));
P_mean  = int_FP/int_F;

% mean melt fraction - Asimow 2001
int_F  = sum(((F(2:end)+F(1:end-1))/2).*abs(P(2:end)-P(1:end-1)));
F_mean = int_F/abs(P(end)-P(1));

% mean major composition
FC = C_fluid_major.*F; FC(1,:) = 0; FC(isnan(FC)) = 0;
int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(P(2:end)-P(1:end-1)));
int_F  = sum(((F(2:end)+F(1:end-1))/2).*abs(P(2:end)-P(1:end-1)));
C_mean_major = int_FC/int_F;

% mean REE composition
FC = C_fluid.*F; FC(1,:) = 0; FC(isnan(FC)) = 0;
int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(P(2:end)-P(1:end-1)));
int_F  = sum(((F(2:end)+F(1:end-1))/2).*abs(P(2:end)-P(1:end-1)));
C_mean_REE = int_FC/int_F;


% mean crustal thickness
g = 9.81;
density_c = 2620;
crust_mean = int_F*1e5/(g*density_c);

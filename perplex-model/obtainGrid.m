function [grid_z] = obtainGrid(M_out,rho_m,g)

P_grid = [M_out(:,1); 0]*1e5;
rho_s  = M_out(end:-1:1,10);
rho_s(1) = rho_m;

DeltaP = P_grid(end-1:-1:1)-P_grid(end:-1:2);
dh = DeltaP./(g*rho_s);

grid_z = flipud(cumsum(dh));
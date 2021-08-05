function [M_out,F_out,grid_z,PTF] = refinement1D(M_out,F_out,ztop,parameters)


index = [1:2:2*size(F_out,1)];
index_interp = [2:2:index(end)-1];

F_out_interp = 0.5*(F_out(1:end-1,:) + F_out(2:end,:));
M_out_interp = 0.5*(M_out(1:end-1,:) + M_out(2:end,:));

F_out_big = zeros(length([index index_interp]),size(F_out,2));
M_out_big = zeros(length([index index_interp]),size(M_out,2));

F_out_big(index,:) = F_out;
F_out_big(index_interp,:) = F_out_interp;
M_out_big(index,:) = M_out;
M_out_big(index_interp,:) = M_out_interp;

F_out = F_out_big;
M_out = M_out_big;

%% second level of refinement for analitical solutions

index = [1:2:2*size(F_out,1)];
index_interp = [2:2:index(end)-1];

F_out_interp = 0.5*(F_out(1:end-1,:) + F_out(2:end,:));
M_out_interp = 0.5*(M_out(1:end-1,:) + M_out(2:end,:));

F_out_big = zeros(length([index index_interp]),size(F_out,2));
M_out_big = zeros(length([index index_interp]),size(M_out,2));

F_out_big(index,:) = F_out;
F_out_big(index_interp,:) = F_out_interp;
M_out_big(index,:) = M_out;
M_out_big(index_interp,:) = M_out_interp;

F_out = F_out_big;
M_out = M_out_big;
%%

P_grid = [M_out(:,1); 0]*1e5;
rho_s  = M_out(end:-1:1,10);
rho_s(1) = parameters.model.rho_m;

g = 9.81;
DeltaP = P_grid(end-1:-1:1)-P_grid(end:-1:2);
dh = DeltaP./(g*rho_s);

grid_z = flipud(cumsum(dh));

% keyboard
[ind_zTop, ~] = find(ztop>grid_z);
if length(ind_zTop)==length(grid_z)
    M_out = [M_out(1,:); M_out];
    F_out = [F_out(1,:); F_out];
    grid_z = [ztop; grid_z];
    [ind_zTop, ~] = find(ztop>grid_z);
end
DeltaZ = grid_z(ind_zTop(1)-1)-grid_z(ind_zTop(1));
weightA = abs(grid_z(ind_zTop(1)-1)-ztop)/DeltaZ;
weightB = abs(grid_z(ind_zTop(1))-ztop)/DeltaZ;
% Interpolate for depth
M_out(ind_zTop(1),:) = M_out(ind_zTop(1)-1,:)*weightA + M_out(ind_zTop(1),:)*weightB;
F_out(ind_zTop(1),:) = F_out(ind_zTop(1)-1,:)*weightA + F_out(ind_zTop(1),:)*weightB;
grid_z(ind_zTop(1),:) = ztop;

% Remove shallower nodes
if length(ind_zTop)>1
    M_out(ind_zTop(2):end,:) = [];
    F_out(ind_zTop(2):end,:) = [];
    grid_z(ind_zTop(2):end) = [];
end

% Remove first nodes (deep nodes) with only solid
ind    = find(M_out(:,12)>0);
if isempty(ind)
    ind = length(M_out(:,12));
    %     keyboard
end
ind_a  = ind(1);
M_out(1:ind_a-2,:)  = [];
F_out(1:ind_a-2,:)  = [];
grid_z(1:ind_a-2,:) = [];


PTF    = [M_out(:,1),M_out(:,2)-273.15,M_out(:,11)];
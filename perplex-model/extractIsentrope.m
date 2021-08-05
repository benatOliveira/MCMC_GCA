function [M_out, F_out, grid_z, PTF] = extractIsentrope(index_Na2O,index_Al2O3,index_Tp,val_Na2O,val_Al2O3,val_Tp,Na2O,Al2O3,Tp,ztop)

index_Na2O_min  = min(index_Na2O);
index_Na2O_max  = max(index_Na2O);
index_Al2O3_min = min(index_Al2O3);
index_Al2O3_max = max(index_Al2O3);
index_Tp_min    = min(index_Tp);
index_Tp_max    = max(index_Tp);
val_Na2O_min  = min(val_Na2O);
val_Na2O_max  = max(val_Na2O);
val_Al2O3_min = min(val_Al2O3);
val_Al2O3_max = max(val_Al2O3);
val_Tp_min    = min(val_Tp);
val_Tp_max    = max(val_Tp);

DeltaNa2O = (Na2O-val_Na2O_min)/(val_Na2O_max-val_Na2O_min);
DeltaAl2O3 = (Al2O3-val_Al2O3_min)/(val_Al2O3_max-val_Al2O3_min); 
DeltaTp = (Tp-val_Tp_min)/(val_Tp_max-val_Tp_min); 

load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_min,index_Tp_min))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_000       = F_out_smooth;
M_out_000       = M_out_smooth;
P_grid_000      = M_out_smooth(:,1);
grid_z_000      = grid_z;
parameters_000  = parameters;
length_000      = size(F_out_000,1);

if DeltaNa2O >= DeltaAl2O3
load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_min,index_Tp_min))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_100       = F_out_smooth;
M_out_100       = M_out_smooth;
P_grid_100      = M_out_smooth(:,1);
grid_z_100      = grid_z;
parameters_100  = parameters;
length_100      = size(F_out_100,1);
else
load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_max,index_Tp_min))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_010       = F_out_smooth;
M_out_010       = M_out_smooth;
P_grid_010      = M_out_smooth(:,1);
grid_z_010      = grid_z;
parameters_010  = parameters;
length_010      = size(F_out_010,1);
end

load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_max,index_Tp_min))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_110       = F_out_smooth;
M_out_110       = M_out_smooth;
P_grid_110      = M_out_smooth(:,1);
grid_z_110      = grid_z;
parameters_110  = parameters;
length_110      = size(F_out_110,1);

load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_min,index_Tp_max))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_001       = F_out_smooth;
M_out_001       = M_out_smooth;
P_grid_001      = M_out_smooth(:,1);
grid_z_001      = grid_z;
parameters_001  = parameters;
length_001      = size(F_out_001,1);

if DeltaNa2O >= DeltaAl2O3
load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_min,index_Tp_max))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_101       = F_out_smooth;
M_out_101       = M_out_smooth;
P_grid_101      = M_out_smooth(:,1);
grid_z_101      = grid_z;
parameters_101  = parameters;
length_101      = size(F_out_101,1);
else
load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_max,index_Tp_max))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_011       = F_out_smooth;
M_out_011       = M_out_smooth;
P_grid_011      = M_out_smooth(:,1);
grid_z_011      = grid_z;
parameters_011  = parameters;
length_011      = size(F_out_011,1);
end

load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_max,index_Tp_max))
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end); 
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
    %% F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    if X(1)<1e-10
        F = 0*F0;
    else
        F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
    end
    F_out_smooth(index,:) = F;
    end
    %%
F_out_111       = F_out_smooth;
M_out_111       = M_out_smooth;
P_grid_111      = M_out_smooth(:,1);
grid_z_111      = grid_z;
parameters_111  = parameters;
length_111      = size(F_out_111,1);

if DeltaNa2O >= DeltaAl2O3
    max_length = max([length_000 length_100 length_110 length_001 length_101 length_111]) ;
    min_length = min([length_000 length_100 length_110 length_001 length_101 length_111]) ;
    if max_length~=length_000; F_out_000 = [repmat(F_out_000(1,:),max_length-length_000,1); F_out_000]; end
    if max_length~=length_100; F_out_100 = [repmat(F_out_100(1,:),max_length-length_100,1); F_out_100]; end
    if max_length~=length_110; F_out_110 = [repmat(F_out_110(1,:),max_length-length_110,1); F_out_110]; end
    if max_length~=length_001; F_out_001 = [repmat(F_out_001(1,:),max_length-length_001,1); F_out_001]; end
    if max_length~=length_101; F_out_101 = [repmat(F_out_101(1,:),max_length-length_101,1); F_out_101]; end
    if max_length~=length_111; F_out_111 = [repmat(F_out_111(1,:),max_length-length_111,1); F_out_111]; end
    if max_length~=length_000; M_out_000 = [repmat(M_out_000(1,:),max_length-length_000,1); M_out_000]; end
    if max_length~=length_100; M_out_100 = [repmat(M_out_100(1,:),max_length-length_100,1); M_out_100]; end
    if max_length~=length_110; M_out_110 = [repmat(M_out_110(1,:),max_length-length_110,1); M_out_110]; end
    if max_length~=length_001; M_out_001 = [repmat(M_out_001(1,:),max_length-length_001,1); M_out_001]; end
    if max_length~=length_101; M_out_101 = [repmat(M_out_101(1,:),max_length-length_101,1); M_out_101]; end
    if max_length~=length_111; M_out_111 = [repmat(M_out_111(1,:),max_length-length_111,1); M_out_111]; end 
else
    max_length = max([length_000 length_010 length_110 length_001 length_011 length_111]) ;
    min_length = min([length_000 length_010 length_110 length_001 length_011 length_111]) ;
    if max_length~=length_000; F_out_000 = [repmat(F_out_000(1,:),max_length-length_000,1); F_out_000]; end
    if max_length~=length_010; F_out_010 = [repmat(F_out_010(1,:),max_length-length_010,1); F_out_010]; end
    if max_length~=length_110; F_out_110 = [repmat(F_out_110(1,:),max_length-length_110,1); F_out_110]; end
    if max_length~=length_001; F_out_001 = [repmat(F_out_001(1,:),max_length-length_001,1); F_out_001]; end
    if max_length~=length_011; F_out_011 = [repmat(F_out_011(1,:),max_length-length_011,1); F_out_011]; end
    if max_length~=length_111; F_out_111 = [repmat(F_out_111(1,:),max_length-length_111,1); F_out_111]; end
    if max_length~=length_000; M_out_000 = [repmat(M_out_000(1,:),max_length-length_000,1); M_out_000]; end
    if max_length~=length_010; M_out_010 = [repmat(M_out_010(1,:),max_length-length_010,1); M_out_010]; end
    if max_length~=length_110; M_out_110 = [repmat(M_out_110(1,:),max_length-length_110,1); M_out_110]; end
    if max_length~=length_001; M_out_001 = [repmat(M_out_001(1,:),max_length-length_001,1); M_out_001]; end
    if max_length~=length_011; M_out_011 = [repmat(M_out_011(1,:),max_length-length_011,1); M_out_011]; end
    if max_length~=length_111; M_out_111 = [repmat(M_out_111(1,:),max_length-length_111,1); M_out_111]; end
end

% if max_length ~= min_length
%     keyboard
% end

if DeltaNa2O >= DeltaAl2O3
    F_out = F_out_000 + (F_out_100-F_out_000)*DeltaNa2O + (F_out_110-F_out_100)*DeltaAl2O3 + (F_out_001-F_out_000)*DeltaTp + ...
                        (F_out_101-F_out_001-F_out_100+F_out_000)*DeltaNa2O*DeltaTp + (F_out_111-F_out_101-F_out_110+F_out_100)*DeltaAl2O3*DeltaTp;
    
    M_out = M_out_000 + (M_out_100-M_out_000)*DeltaNa2O + (M_out_110-M_out_100)*DeltaAl2O3 + (M_out_001-M_out_000)*DeltaTp + ...
                        (M_out_101-M_out_001-M_out_100+M_out_000)*DeltaNa2O*DeltaTp + (M_out_111-M_out_101-M_out_110+M_out_100)*DeltaAl2O3*DeltaTp;
else
    F_out = F_out_000 + (F_out_110-F_out_010)*DeltaNa2O + (F_out_010-F_out_000)*DeltaAl2O3 + (F_out_001-F_out_000)*DeltaTp + ...
                        (F_out_111-F_out_011-F_out_110+F_out_010)*DeltaNa2O*DeltaTp + (F_out_011-F_out_001-F_out_010+F_out_000)*DeltaAl2O3*DeltaTp;
    
    M_out = M_out_000 + (M_out_110-M_out_010)*DeltaNa2O + (M_out_010-M_out_000)*DeltaAl2O3 + (M_out_001-M_out_000)*DeltaTp + ...
                        (M_out_111-M_out_011-M_out_110+M_out_010)*DeltaNa2O*DeltaTp + (M_out_011-M_out_001-M_out_010+M_out_000)*DeltaAl2O3*DeltaTp;
end
 

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


    

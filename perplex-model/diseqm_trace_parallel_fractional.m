function [C_solid, C_fluid, C_extracted, C_ave, C_solid_major, C_fluid_major, C_extracted_major, C_ave_major, C_Ol, C_Cpx, C_Opx, C_Grt, C_Sp, C_Pl, w_j, phi_j, v_s, v_f, trace_fail] = diseqm_trace_parallel_fractional(M_out,F_out,grid_z,parameters,C_source_REE)
% DISEQM_TRACE_PARALLEL_FRACTIONAL
%
% [C_solid, C_fluid, C_extracted, C_ave, C_solid_major, C_fluid_major, C_extracted_major, C_ave_major, C_Ol, C_Cpx, C_Opx, C_Grt, C_Sp, C_Pl, w_j, phi_j] = diseqm_trace_parallel_fractional(M_out,F_out,grid_z,parameters)
%
% INPUT:
%
% M_out
% F_out
% grid_z
% parameters
%
% OUTPUT:
%
% C_solid:          Trace element composition within the solid (normalized)
% C_fluid:          Trace element composition within the liquid allowed to equilibrate  (normalized)
% C_extracted:      Trace element composition within the extracted liquid (normalized)
% C_ave:            Trace element composition within the total liquid in average (normalized)
%       Columns:    1    2    3    4    5    6    7    8    9    10   11   12   13
%                   La - Ce - Pr - Nd - Sm - Eu - Gd - Tb - Dy - Ho - Er - Yb - Lu

% C_solid_major:    Major element composition within the solid
% C_fluid_major:    Major element composition within the liquid allowed to equilibrate
% C_extracted_major Major element composition within the extracted liquid
% C_ave_major:      Major element composition within the total liquid in average (normalized)
%       Columns:    1      2       3     4     5     6     7      8       9
%                   SiO2 - Al2O3 - FeO - MnO - MgO - CaO - Na2O - Cr2O3 - TiO2

% w_j:              weight percentage of each phase [thermodynamic phase weight]/ [total (solid+liquid) weight]
% phi_j:            volume percentage of each phase [thermodynamic phase volume]/ [total (solid+liquid) volume]

%% set-up

nz = size(M_out); nz = nz(1);                       % number of intervals/nodes

P = M_out(:,1);
T = M_out(:,2);
F = F_out(:,1);
H = grid_z(1)-grid_z(end);
U = 2*((grid_z-grid_z(end))/H) - ((grid_z-grid_z(end))/H).^2;
if length(grid_z)==2
    U=ones(size(U,1),1); % velocity weight => 1
end
W0 = parameters.model.W0/(60*60*24*365)/100;           % solid upwelling velocity [m/s]

phi_inst = M_out(:,12)./100;
fT = parameters.model.fT;                               % residual porosity - traces
g = parameters.model.g;                               % gravitational acceleration [m/s2]
k0 = parameters.model.k0;                                         % permeability 1e-6 to 1e-11
mu = parameters.model.mu_l;                                           % fluid viscosity [Pa s]
phi = zeros(nz,1);
phi_eq = zeros(nz,1);
phi_col = 0*ones(nz,1);
mu = mu*ones(nz,1);
k0 = k0*ones(nz,1);
n = parameters.model.n;
rho_s = M_out(:,10);                                 % bulk residue
rho_l = parameters.model.rho_l;
rho_l = rho_l*ones(nz,1);
v_s = W0*ones(nz,1);
v_f = W0*ones(nz,1);
S_int = 0*ones(nz,1);
S_dis = 0*ones(nz,1);
Drho = rho_s-rho_l;
beta = parameters.model.beta;                        % geometry of crystals: 5 = sphere, 4 = cylinder, 1 = plane sheet
trace_fail = 0;


%% Partition Coefficients, Diffusion

traceImport = importTrace('trace_elements_properties.txt');
[trace]     = createTrace(traceImport);
% REE mineral partition coefficients -- cpx, Gibson&Geist (2010)
%  1     2     3     4     5     6     7     8     9     10    11    12    13
%  La    Ce    Pr    Nd    Sm    Eu    Gd    Tb    Dy    Ho    Er    Yb    Lu

% olivine, cpx, opx, garnet, spinel
K_j = cell(nz,13);

% olivine, cpx, opx, garnet, spinel, pl
K_j(1:nz,1) = {trace.K(1,1:end)};
K_j(1:nz,2) = {trace.K(2,1:end)};
K_j(1:nz,3) = {trace.K(3,1:end)};
K_j(1:nz,4) = {trace.K(4,1:end)};
K_j(1:nz,5) = {trace.K(5,1:end)};
K_j(1:nz,6) = {trace.K(6,1:end)};
K_j(1:nz,7) = {trace.K(7,1:end)};
K_j(1:nz,8) = {trace.K(8,1:end)};
K_j(1:nz,9) = {trace.K(9,1:end)};
K_j(1:nz,10) = {trace.K(10,1:end)};
K_j(1:nz,11) = {trace.K(11,1:end)};
K_j(1:nz,12) = {trace.K(12,1:end)};
K_j(1:nz,13) = {trace.K(13,1:end)};

% check for P-T-C dependent partition coefficients
if parameters.model.part == 1
    K_j = changePartitionCoefficient(K_j,M_out(:,48:56),M_out(:,57:65),M_out(:,66:74),P,T,M_out(:,12:18),grid_z);
end


%% Diffusion coefficient
R = 8.314;                      % [J/(mol K)]; gas constant

% % olivine, cpx, opx, garnet, spinel, pl
epsilon = trace.E.*1e3;
V = trace.V;
factor_La = exp((-epsilon(1,1:end) + P*1e-5*V(1,1:end))./(R*(T)));
factor_Ce = exp((-epsilon(2,1:end) + P*1e-5*V(2,1:end))./(R*(T)));
factor_Pr = exp((-epsilon(3,1:end) + P*1e-5*V(3,1:end))./(R*(T)));
factor_Nd = exp((-epsilon(4,1:end) + P*1e-5*V(4,1:end))./(R*(T)));
factor_Sm = exp((-epsilon(5,1:end) + P*1e-5*V(5,1:end))./(R*(T)));
factor_Eu = exp((-epsilon(6,1:end) + P*1e-5*V(6,1:end))./(R*(T)));
factor_Gd = exp((-epsilon(7,1:end) + P*1e-5*V(7,1:end))./(R*(T)));
factor_Tb = exp((-epsilon(8,1:end) + P*1e-5*V(8,1:end))./(R*(T)));
factor_Dy = exp((-epsilon(9,1:end) + P*1e-5*V(9,1:end))./(R*(T)));
factor_Ho = exp((-epsilon(10,1:end) + P*1e-5*V(10,1:end))./(R*(T)));
factor_Er = exp((-epsilon(11,1:end) + P*1e-5*V(11,1:end))./(R*(T)));
factor_Yb = exp((-epsilon(12,1:end) + P*1e-5*V(12,1:end))./(R*(T)));
factor_Lu = exp((-epsilon(13,1:end) + P*1e-5*V(13,1:end))./(R*(T)));

D_La = trace.D(1,1:end).*factor_La;
D_Ce = trace.D(2,1:end).*factor_Ce;
D_Pr = trace.D(3,1:end).*factor_Pr;
D_Nd = trace.D(4,1:end).*factor_Nd;
D_Sm = trace.D(5,1:end).*factor_Sm;
D_Eu = trace.D(6,1:end).*factor_Eu;
D_Gd = trace.D(7,1:end).*factor_Gd;
D_Tb = trace.D(8,1:end).*factor_Tb;
D_Dy = trace.D(9,1:end).*factor_Dy;
D_Ho = trace.D(10,1:end).*factor_Ho;
D_Er = trace.D(11,1:end).*factor_Er;
D_Yb = trace.D(12,1:end).*factor_Yb;
D_Lu = trace.D(13,1:end).*factor_Lu;

Dj = cell(nz,13);
Dj(1:nz,1)  = mat2cell(D_La,ones(1,nz),6);
Dj(1:nz,2)  = mat2cell(D_Ce,ones(1,nz),6);
Dj(1:nz,3)  = mat2cell(D_Pr,ones(1,nz),6);
Dj(1:nz,4)  = mat2cell(D_Nd,ones(1,nz),6);
Dj(1:nz,5)  = mat2cell(D_Sm,ones(1,nz),6);
Dj(1:nz,6)  = mat2cell(D_Eu,ones(1,nz),6);
Dj(1:nz,7)  = mat2cell(D_Gd,ones(1,nz),6);
Dj(1:nz,8)  = mat2cell(D_Tb,ones(1,nz),6);
Dj(1:nz,9)  = mat2cell(D_Dy,ones(1,nz),6);
Dj(1:nz,10) = mat2cell(D_Ho,ones(1,nz),6);
Dj(1:nz,11) = mat2cell(D_Er,ones(1,nz),6);
Dj(1:nz,12) = mat2cell(D_Yb,ones(1,nz),6);
Dj(1:nz,13) = mat2cell(D_Lu,ones(1,nz),6);

dj = cell(nz,13);
dj(:,:) = {trace.radii(1:end)};                             % typical/average grain size for each mineral phase [m]

% Define parameter for diffusion
R_j = (3*beta.*cell2mat(Dj))./(cell2mat(dj).^2);  % parameter for diffusion and grain size
R_j(R_j>1e-11) = 1e-11;
if parameters.model.diff ~= 1 % check if constant R required
    R_j = (1e-12)*ones(size(R_j));
    
end
R_j = mat2cell(R_j,ones(1,nz),ones(1,13)*6);


%% Run REE model
if isnan(sum(phi_inst)) || sum(phi_inst)==0
    La = 0; Ce = 0; Pr = 0; Nd = 0; Sm = 0; Eu = 0; Gd = 0; Tb = 0; Dy = 0; Ho = 0; Er = 0; Yb = 0; Lu = 0;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for phi, v_s and v_f
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ind = find(phi_inst>0); ind_a = ind(1); ind_o = ind(end);
    grid_dist = grid_z(1:end-1)-grid_z(2:end);
    rho0_s = rho_s(ind_a);
    phi_crit = fT;        % critical porosity (retained melt in matrix allowed to equilibrate)
    
    for index=ind_a:nz
        if n == 3
            p = [1 ,...
                -2 ,...
                1 ,...
                0 , ...
                (-F(index)*W0*mu(index)*rho_l(index)*rho0_s + F(index)*W0*mu(index)*rho0_s*rho_s(index) + W0*mu(index)*rho_l(index)*rho0_s)/(Drho(index)*g*k0(index)*rho_l(index)*rho_s(index)) , ...
                (- F(index)*W0*mu(index)*rho0_s*rho_s(index))/(Drho(index)*g*k0(index)*rho_l(index)*rho_s(index))];
            sol = roots(p);
            phi_aux = sol(imag(sol)==0);
        elseif n == 2
            p = [1 ,...
                -2 ,...
                1 ,...
                (-F(index)*W0*mu(index)*rho_l(index)*rho0_s + F(index)*W0*mu(index)*rho0_s*rho_s(index) + W0*mu(index)*rho_l(index)*rho0_s)/(Drho(index)*g*k0(index)*rho_l(index)*rho_s(index)) , ...
                (- F(index)*W0*mu(index)*rho0_s*rho_s(index))/(Drho(index)*g*k0(index)*rho_l(index)*rho_s(index))];
            sol = roots(p);
            phi_aux = sol(imag(sol)==0);
        else
            error('permeability exponent not implemented')
            
        end
        
        phi(index) = phi_aux;
        if phi(index)>phi_crit && phi_crit~=0   %  case where there is some melt fraction to suck
            phi_col(index) = phi(index)-phi_crit;
            phi_eq(index)= phi_crit;
            S_int(index) = (k0(index)*phi_crit^n*rho_l(index)*(Drho(index)*g*(phi_col(index) + phi_crit - 1) + (mu(index)*((F(index)*W0*rho0_s)/(phi_crit*rho_l(index)) - (W0*rho0_s*(F(index) - 1))/(rho_s(index)*(phi_col(index) + phi_crit - 1))))/(k0(index)*phi_crit^(n-1))))/mu(index);
            S_dis(index) = (S_int(index)-S_int(index-1))/grid_dist(index-1);
        elseif phi_crit == 0    %  case where all the melt fraction needs to be sucked
            phi_col(index) = phi(index)-phi_crit;
            phi_eq(index)= phi_crit;
            S_int(index) = F(index)*W0*rho0_s;
            S_dis(index) = (S_int(index)-S_int(index-1))/grid_dist(index-1);
        else                    %  case where no melt fraction needs to be sucked
            phi_eq(index)= phi(index);
            S_int(index) = 0;
            S_dis(index) = (S_int(index)-S_int(index-1))/grid_dist(index-1);
        end
        
        v_s(index) = ((1-F(index))*rho0_s*W0)/((1-phi_eq(index)-phi_col(index))*rho_s(index));
        v_f(index) = (F(index)*rho0_s*W0)/((phi_eq(index)+phi_col(index))*rho_l(index));
        if phi(index)==0
            v_f(index) = v_s(index);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties for the ODE solvers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind = find(phi>0); ind_a = ind(1); ind_o = ind(end);
    
    phi_s = 1 - phi;
    rho_t = phi.*rho_l + phi_s.*rho_s;
    w_s = phi_s.*rho_s./rho_t;
    w_f = 1- w_s;
    tol_phi = 1e-3; % per 1
    
    TP_weight   = M_out(:,13:18);
    TP_volume   = M_out(:,20:25);
    
    TP_weight(TP_weight<tol_phi | TP_volume<tol_phi) = 0;
    TP_volume(TP_weight<tol_phi | TP_volume<tol_phi) = 0;
    
    % olivine, cpx, opx, garnet, spinel, pl
    TP_weight   = TP_weight./repmat(sum(TP_weight,2),1,6)*100;
    TP_volume   = TP_volume./repmat(sum(TP_volume,2),1,6)*100;
    phi_j       = TP_volume(:,[1 2 3 4 5 6])./100.*repmat(1-phi,1,6);         % [thermodynamic phase volume]/ [total (solid+melt) volume]
    phi_j_solid = TP_volume(:,[1 2 3 4 5 6])./100;                            % [thermodynamic phase volume]/ [solid volume]
    w_j         = TP_weight(:,[1 2 3 4 5 6])./100.*repmat(w_s,1,6);           % [thermodynamic phase weight]/ [total (solid+melt) weight]
    w_j_solid   = TP_weight(:,[1 2 3 4 5 6])./100;                            % [thermodynamic phase weight]/ [solid weight]
    rho_j       = M_out(:,27:32);
    
    % plot some intermidiate results - uncomment if needed
    %     figure(5);clf(5);
    %     subplot(1,5,1)
    %     plot(100*phi_j,M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('\phi_j (%wt)');legend('Ol','Cpx','Opx','Gt','Sp','Pl')
    %     subplot(1,5,2)
    %     plot(v_f*100*(365*24*60*60),M_out(:,1)*-1e-3);hold on;plot(v_s*100*(365*24*60*60),M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('velocity (cm/yr)');legend('liquid','solid')
    %     subplot(1,5,3)
    %     plot(M_out(:,2) - 273.15,M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('T (C)');
    %     subplot(1,5,4)
    %     plot(100*S_int/(W0*rho0_s),M_out(:,1)*-1e-3);hold on;plot(100*F,M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('\int S^{modi} and F (%wt)');legend('\int S^*','F')
    %     subplot(1,5,5)
    %     scatter(100*phi,M_out(:,1)*-1e-3,'filled');hold on;plot(100*phi_col,M_out(:,1)*-1e-3);plot(100*phi_eq,M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('\phi (%wt)');legend('\phi','\phi_{col}','\phi_{eq}')
    
    
    % Compute Gamma out of F
    F_out(isnan(F_out)) = 0;  F_out = F_out - F_out(ind_a-1,:);
    Gamma = -rho0_s*W0.*(F_out(1:end-1,:)-F_out(2:end,:))./repmat(grid_dist,1,size(F_out,2));    % negative sign is to take into account the z axis
    Gamma = [0 0 0 0 0 0 0; Gamma];
    
    % Rearraange Gamma for ODE
    Gamma_ode = -Gamma(:,[2 3 4 5 6 7]);   % positive ==> melting phase
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve MELTS MAJOR OXIDE CHEMISTRY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initial conditions
    %       Columns:    1      2       3     4     5     6     7      8       9
    %                   SiO2 - Al2O3 - FeO - MnO - MgO - CaO - Na2O - Cr2O3 - TiO2
    C_solid_aux   = M_out(:,39:47)./repmat(sum(M_out(:,39:47),2),1,9); C_solid = C_solid_aux;
    C_melt_aux    = M_out(:,48:56)./repmat(sum(M_out(:,48:56),2),1,9); C_melt  = C_melt_aux;
    
    C_fluid_major        = 0*[C_solid];
    C_extracted_major    = 0*[C_solid];
    
    
    if ind_a ~= ind_o
        
        % interpolate values for smooth parameters
        C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
        C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:); C_solid(C_solid<0) = 0; C_solid = C_solid./repmat(sum(C_solid,2),1,9);
        C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
        C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:); C_melt(C_melt<0) = 0; C_melt = C_melt./repmat(sum(C_melt,2),1,9);
        
        % Numerical solution for extracted melt composition - Eq. 2.2 Zou
        C_melt_zou = C_melt;
        C_melt_zou(isnan(C_melt_zou)) = 0;
        int_C = cumsum(((C_melt_zou(2:end,:)+C_melt_zou(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)));
        int_C = cumsum(((C_melt_zou(2:end,:)+C_melt_zou(1:end-1,:))/2).*abs(F(2:end)-F(1:end-1)));
        F_half = (F(2:end)+F(1:end-1))/2;
        OX_extr = int_C./F_half; OX_extr = OX_extr./(sum(OX_extr,2));
        OX_extr(isnan(OX_extr)) = 0; OX_extr = [0*OX_extr(1,:); OX_extr];
        
        C = OX_extr;
        
    else
        
        C = [C_melt(ind_a,:);...
            C_melt(ind_a,:)];
    end
    
    C_fluid_major(ind_a:ind_o,:) = C_melt(ind_a:ind_o,:);
    C_solid_major = C_solid;
    C_extracted_major = C;
    
    C_total_fluid_major = (repmat(phi_eq,1,size(C_fluid_major,2)).*C_fluid_major + repmat(phi_col,1,size(C_extracted_major,2)).*C_extracted_major)./repmat(phi_eq+phi_col,1,size(C_extracted_major,2));
    
    % Numerical solution for averaged melt composition - e.g. Ito 2005
    % Eq. 1
    % weighted with F
    FC = C_total_fluid_major.*F.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((F(2:end).*U(2:end)+F(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_major_F = int_FC/int_F;
    % weighted with \phi
    FC = C_total_fluid_major.*phi.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((phi(2:end).*U(2:end)+phi(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_major_phi = int_FC/int_F;
    % weighted with \phi v_f
    FC = C_total_fluid_major.*phi.*v_f.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((phi(2:end).*v_f(2:end).*U(2:end)+phi(1:end-1).*v_f(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_major_phi_v = int_FC/int_F;
    
    C_ave_major = C_ave_major_F;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVE TRACE ELEMENT CHEMISTRY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%% when we invert for source composition%%%%%%%%%%%%%%%%%
    %  1       2       3       4       5       6       7       8       9       10      11      12      13
    %  La      Ce      Pr      Nd      Sm      Eu      Gd      Tb      Dy      Ho      Er      Yb      Lu
    C_source=C_source_REE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Steady-state ODE
    for index_par = 1:length(C_source)
        %% Trace
        Kj_element = cell2mat(K_j(:,index_par));
        Rj_element = cell2mat(R_j(:,index_par));
        
        Cf_init = repmat(C_source(index_par),length(grid_z),1)./(sum(w_j.*Kj_element,2)+w_f);
        Cs_init = repmat(Cf_init,1,size(Kj_element,2)).*Kj_element; Cs_init(w_j==0) = 0;
        La        = [0*Cf_init Cs_init 0*Cf_init];  % initialize
        
        init_cond = [Cf_init(ind_a),Cs_init(ind_a,:),Cf_init(ind_a)]';                     % [fluid,mineral1,mineral2,mineral3,...,mineralN,extractedFluid]
        if ind_a ~= ind_o
            step_nz = -[0 ;grid_z(2:end)-grid_z(1:end-1)];
            XX = zeros(length(grid_z),length(init_cond));
            XX(1:ind_a,:) = repmat(init_cond',ind_a,1);
            for index_nz= ind_a+1:ind_o
                if phi_col(index_nz) == 0 
                    A11 = sum(Gamma_ode(index_nz,:)) + rho_s(index_nz)*(1-phi_col(index_nz)-phi_eq(index_nz))*sum(w_j(index_nz,:).*Rj_element(index_nz,:).*Kj_element(index_nz,:));
                    A12 = - Gamma_ode(index_nz,:) - rho_s(index_nz)*(1-phi_col(index_nz)-phi_eq(index_nz))*w_j(index_nz,:).*Rj_element(index_nz,:);
                    A13 = 0;
                    A21 = - step_nz(index_nz)/(2*v_s(index_nz))*(Rj_element(index_nz,:)'.*Kj_element(index_nz,:)');
                    A22 = diag(ones(6,1)) + diag(step_nz(index_nz)/(2*v_s(index_nz))*Rj_element(index_nz,:));
                    A23 = zeros(6,1);
                    A31 = - step_nz(index_nz)/(2*rho_l(index_nz)*phi_col(index_nz)*v_f(index_nz))*S_dis(index_nz);
                    A32 = zeros(1,6);
                    A33 = 1 + step_nz(index_nz)/(2*rho_l(index_nz)*phi_col(index_nz)*v_f(index_nz))*S_dis(index_nz);
                    
                    b1 = 0;
                    b2 = XX(index_nz-1,2:end-1)' - step_nz(index_nz)/(2*v_s(index_nz-1))*Rj_element(index_nz-1,:)'.*XX(index_nz-1,2:end-1)' + ...
                                                 + step_nz(index_nz)/(2*v_s(index_nz-1))*Rj_element(index_nz-1,:)'.*Kj_element(index_nz-1,:)'.*XX(index_nz-1,1)';
                    b3 = XX(index_nz-1,end)      - step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_col(index_nz-1)*v_f(index_nz-1))*S_dis(index_nz-1).*XX(index_nz-1,end) + ...
                                                 + step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_col(index_nz-1)*v_f(index_nz-1))*S_dis(index_nz-1).*XX(index_nz-1,1);
                    
                elseif phi_eq(index_nz) == 0
                    A11 = sum(Gamma_ode(index_nz,:)) + rho_s(index_nz)*(1-phi_col(index_nz)-phi_eq(index_nz))*sum(w_j(index_nz,:).*Rj_element(index_nz,:).*Kj_element(index_nz,:));
                    A12 = - Gamma_ode(index_nz,:) - rho_s(index_nz)*(1-phi_col(index_nz)-phi_eq(index_nz))*w_j(index_nz,:).*Rj_element(index_nz,:);
                    A13 = 0;
                    A21 = - step_nz(index_nz)/(2*v_s(index_nz))*(Rj_element(index_nz,:)'.*Kj_element(index_nz,:)');
                    A22 = diag(ones(6,1)) + diag(step_nz(index_nz)/(2*v_s(index_nz))*Rj_element(index_nz,:));
                    A23 = zeros(6,1);
                    A31 = - step_nz(index_nz)/(2*rho_l(index_nz)*phi_col(index_nz)*v_f(index_nz))*S_dis(index_nz);
                    A32 = zeros(1,6);
                    A33 = 1 + step_nz(index_nz)/(2*rho_l(index_nz)*phi_col(index_nz)*v_f(index_nz))*S_dis(index_nz);
                    
                    b1 =                         - rho_s(index_nz-1)*(1-phi_col(index_nz-1)-phi_eq(index_nz-1))*sum(w_j(index_nz-1,:).*Rj_element(index_nz-1,:).*Kj_element(index_nz-1,:)).*XX(index_nz-1,1) + ...
                        - sum(Gamma_ode(index_nz-1,:)).*XX(index_nz-1,1) + ...
                        + rho_s(index_nz-1)*(1-phi_col(index_nz-1)-phi_eq(index_nz-1))*sum(w_j(index_nz-1,:).*Rj_element(index_nz-1,:).*XX(index_nz-1,2:end-1)) + ...
                        + sum(Gamma_ode(index_nz-1,:).*XX(index_nz-1,2:end-1));
                    b1 = 0;
                    b2 = XX(index_nz-1,2:end-1)' - step_nz(index_nz)/(2*v_s(index_nz-1))*Rj_element(index_nz-1,:)'.*XX(index_nz-1,2:end-1)' + ...
                        + step_nz(index_nz)/(2*v_s(index_nz-1))*Rj_element(index_nz-1,:)'.*Kj_element(index_nz-1,:)'.*XX(index_nz-1,1)';
                    b3 = XX(index_nz-1,end)      - step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_col(index_nz-1)*v_f(index_nz-1))*S_dis(index_nz-1).*XX(index_nz-1,end) + ...
                        + step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_col(index_nz-1)*v_f(index_nz-1))*S_dis(index_nz-1).*XX(index_nz-1,1);  b3(S_dis(index_nz-1)==0) = 0;
                    
                else
                    A11 = 1 + step_nz(index_nz)/(2*rho_l(index_nz)*phi_eq(index_nz)*v_f(index_nz))*sum(Gamma_ode(index_nz,:)) + step_nz(index_nz)/(2*rho_l(index_nz)*phi_eq(index_nz)*v_f(index_nz))*rho_s(index_nz)*(1-phi_col(index_nz)-phi_eq(index_nz))*sum(w_j(index_nz,:).*Rj_element(index_nz,:).*Kj_element(index_nz,:));
                    A12 = - step_nz(index_nz)/(2*rho_l(index_nz)*phi_eq(index_nz)*v_f(index_nz))*Gamma_ode(index_nz,:) - step_nz(index_nz)/(2*rho_l(index_nz)*phi_eq(index_nz)*v_f(index_nz))*rho_s(index_nz)*(1-phi_col(index_nz)-phi_eq(index_nz))*w_j(index_nz,:).*Rj_element(index_nz,:);
                    A13 = 0;
                    A21 = - step_nz(index_nz)/(2*v_s(index_nz))*(Rj_element(index_nz,:)'.*Kj_element(index_nz,:)');
                    A22 = diag(ones(6,1)) + diag(step_nz(index_nz)/(2*v_s(index_nz))*Rj_element(index_nz,:));
                    A23 = zeros(6,1);
                    A31 = - step_nz(index_nz)/(2*rho_l(index_nz)*phi_col(index_nz)*v_f(index_nz))*S_dis(index_nz);
                    A32 = zeros(1,6);
                    A33 = 1 + step_nz(index_nz)/(2*rho_l(index_nz)*phi_col(index_nz)*v_f(index_nz))*S_dis(index_nz);
                    
                    
                    
                    b1 = XX(index_nz-1,1)        - step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_eq(index_nz-1)*v_f(index_nz-1))*rho_s(index_nz-1)*(1-phi_col(index_nz-1)-phi_eq(index_nz-1))*sum(w_j(index_nz-1,:).*Rj_element(index_nz-1,:).*Kj_element(index_nz-1,:)).*XX(index_nz-1,1) + ...
                        - step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_eq(index_nz-1)*v_f(index_nz-1))*sum(Gamma_ode(index_nz-1,:)).*XX(index_nz-1,1) + ...
                        + step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_eq(index_nz-1)*v_f(index_nz-1))*rho_s(index_nz-1)*(1-phi_col(index_nz-1)-phi_eq(index_nz-1))*sum(w_j(index_nz-1,:).*Rj_element(index_nz-1,:).*XX(index_nz-1,2:end-1)) + ...
                        + step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_eq(index_nz-1)*v_f(index_nz-1))*sum(Gamma_ode(index_nz-1,:).*XX(index_nz-1,2:end-1));
                    b2 = XX(index_nz-1,2:end-1)' - step_nz(index_nz)/(2*v_s(index_nz-1))*Rj_element(index_nz-1,:)'.*XX(index_nz-1,2:end-1)' + ...
                        + step_nz(index_nz)/(2*v_s(index_nz-1))*Rj_element(index_nz-1,:)'.*Kj_element(index_nz-1,:)'.*XX(index_nz-1,1)';
                    b3 = XX(index_nz-1,end)      - step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_col(index_nz-1)*v_f(index_nz-1))*S_dis(index_nz-1).*XX(index_nz-1,end) + ...
                        + step_nz(index_nz)/(2*rho_l(index_nz-1)*phi_col(index_nz-1)*v_f(index_nz-1))*S_dis(index_nz-1).*XX(index_nz-1,1);  b3(S_dis(index_nz-1)==0) = 0;
                end
                
                %         % scale
                %         A11 = step_nz(index_nz)/(2*v_s(index_nz))*A11;
                %         A12 = step_nz(index_nz)/(2*v_s(index_nz))*A12;
                %         b1 = step_nz(index_nz)/(2*v_s(index_nz))*b1;
                
                ind_remove = [false w_j(index_nz,:)==0 phi_col(index_nz,:)==0];
                ind_comput = ind_remove==0;
                A = [   A11 A12 A13;
                    A21 A22 A23;
                    A31 A32 A33];
                
                b = [   b1; b2; b3 ];
                
                A(:,ind_remove) = [];
                A(ind_remove,:) = [];
                b(ind_remove)   = [];
                
                result = A\b;
                
                if any(isnan(result))
                    trace_fail = 1; 
                end
                
                
                if any(result<0); [result,~,~,flag] = lsqnonneg(A,b); if flag==0; trace_fail = 1; end; end
                
                XX(index_nz,:) = XX(index_nz-1,:);
                XX(index_nz,:) = 0;
                XX(index_nz,ind_comput) = result;
                
            end
        else
            XX = [init_cond';...
                init_cond'];
        end
        
        % Save variables
        La(ind_a:ind_o,:) = XX(ind_a:ind_o,:);
        La_solid = sum(La(:,2:end-1).*w_j_solid,2);
        La_extracted = La(:,end);
        
        REE_inst{index_par}  = La(:,1);
        REE_extr{index_par}  = La_extracted;
        REE_solid{index_par} = La_solid;
        
        % olivine, cpx, opx, garnet, spinel, pl
        REE_ol{index_par}    = La(:,2);
        REE_cpx{index_par}   = La(:,3);
        REE_opx{index_par}   = La(:,4);
        REE_grt{index_par}   = La(:,5);
        REE_sp{index_par}    = La(:,6);
        REE_pl{index_par}    = La(:,7);
        
        
    end
    
    % Normalize
    %                  La      Ce      Pr      Nd      Sm      Eu      Gd      Tb    Dy      Ho      Er      Yb      Lu
    REE_factor =      [237     613     92.8    457     148     56.3    199     36.1  246     54.6    160     161     24.6]/1000; % CI from McDonough & Sun (1995)
    
    C_solid = cell2mat(REE_solid)./repmat(REE_factor,length(REE_solid{1}),1);
    C_fluid = cell2mat(REE_inst)./repmat(REE_factor,length(REE_inst{1}),1);   % reservoir 1 (instant. liquid allowed to equilibrate)
    C_extracted = cell2mat(REE_extr)./repmat(REE_factor,length(REE_extr{1}),1); % reservoir 2 (isolated liquid)
    C_total_fluid = (repmat(phi_eq,1,size(C_fluid,2)).*C_fluid + repmat(phi_col,1,size(C_extracted,2)).*C_extracted)./repmat(phi_eq+phi_col,1,size(C_extracted,2));
    
    
    % Numerical solution for averaged melt composition - Ito 2005
    % Eq. 1
    % weighted with F
    FC = C_total_fluid.*F.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((F(2:end).*U(2:end)+F(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_F = int_FC/int_F;
    % weighted with \phi
    FC = C_total_fluid.*phi.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((phi(2:end).*U(2:end)+phi(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_phi = int_FC/int_F;
    % weighted with \phi v_f
    FC = C_total_fluid.*phi.*v_f.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((phi(2:end).*v_f(2:end).*U(2:end)+phi(1:end-1).*v_f(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_phi_v = int_FC/int_F;
    % weighted with F a la McKenzie
    U=ones(size(U,1),1); % velocity weight => 1
    FC = C_total_fluid.*F.*U; FC(1,:) = 0; FC(isnan(FC)) = 0;
    int_FC = sum(((FC(2:end,:)+FC(1:end-1,:))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    int_F = sum(((F(2:end).*U(2:end)+F(1:end-1).*U(1:end-1))/2).*abs(grid_z(2:end)-grid_z(1:end-1)),1);
    C_ave_FM = int_FC/int_F;
    
    % Choose final form for pooled melt composition...
    C_ave = C_ave_F;
    
    % Some checks for exploratory purposes
    C_ave(1)/C_ave(5);
    C_ave_FM(1)/C_ave_FM(5);
    
    % save thermodynamic phase REE composition
    C_Ol    = cell2mat(REE_ol)./repmat(REE_factor,length(REE_ol{1}),1);
    C_Cpx   = cell2mat(REE_cpx)./repmat(REE_factor,length(REE_cpx{1}),1);
    C_Opx   = cell2mat(REE_opx)./repmat(REE_factor,length(REE_opx{1}),1);
    C_Grt   = cell2mat(REE_grt)./repmat(REE_factor,length(REE_grt{1}),1);
    C_Sp    = cell2mat(REE_sp)./repmat(REE_factor,length(REE_sp{1}),1);
    C_Pl    = cell2mat(REE_pl)./repmat(REE_factor,length(REE_pl{1}),1);
    
    %% Check results and flagged if negative, zero or nan
    tf = isreal(C_ave);
    tnan=isnan(C_ave);
    if (any(C_ave <= 0) || tf == 0 || any(tnan == 1))
        trace_fail = 1;         
    end
    
    %% Plot - uncomment if necessary - check if you are using parallel pool
    
    Colormap=jet(size(C_solid,1));
    
    
% % % %         figure(10);clf(10);
% % % %         subplot(1,6,1)
% % % %         h=semilogy(C_solid');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_s^e (ppm)');xlabel('REE');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,2)
% % % %         h=semilogy(C_total_fluid');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_l^e (ppm)');xlabel('REE'); % instantaneo
% % % %                             hold on;  plot(C_ave','.r', 'MarkerSize',15')
% % % %                             hold on;  plot(C_ave_phi_v','.b', 'MarkerSize',20')
% % % %                             hold on;  plot(C_ave_FM','og', 'MarkerSize',6')
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,3)
% % % %         plot(v_f*100*(365*24*60*60),M_out(:,1)*-1e-3,'.');hold on;plot(v_s*100*(365*24*60*60),M_out(:,1)*-1e-3,'.b');ylabel('P (kbar)');xlabel('velocity (cm/yr)');legend('liquid','solid')
% % % %         subplot(1,6,4)
% % % %         plot(M_out(:,2) - 273.15,M_out(:,1)*-1e-3,'.');ylabel('P (kbar)');xlabel('T (C)');
% % % %         subplot(1,6,5)
% % % %         plot(100*phi_j_solid,M_out(:,1)*-1e-3,'.');ylabel('P (kbar)');xlabel('\phi_j (%wt)');legend('Ol','Cpx','Opx','Gt','Sp')
% % % %         subplot(1,6,6)
% % % %         plot(100*phi_eq,M_out(:,1)*-1e-3,'.');hold on;plot(100*phi_col,M_out(:,1)*-1e-3,'.');plot(100*F,M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('\phi_{eq}, \phi_{col} and F (%wt)');legend('\phi_{eq}','\phi_{col}','F')
% % % %     
% % % % 
% % % %     
% % % %     
% % % %         figure(20);clf(20);
% % % %         subplot(1,6,1)
% % % %         h=semilogy(C_Ol');xlim([1 13]);ylim([1e-6 1e1]);ylabel('C_s^e (ppm)');xlabel('REE - Ol');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,2)
% % % %         h=semilogy(C_Cpx');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_s^e (ppm)');xlabel('REE - Cpx');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,3)
% % % %         h=semilogy(C_Opx');xlim([1 13]);ylim([1e-3 1e1]);ylabel('C_s^e (ppm)');xlabel('REE - Opx');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,4)
% % % %         h=semilogy(C_Sp');xlim([1 13]);ylim([1e-3 1e1]);ylabel('C_s^e (ppm)');xlabel('REE - Sp');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,5)
% % % %         h=semilogy(C_Grt');xlim([1 13]);ylim([1e-6 1e1]);ylabel('C_s^e (ppm)');xlabel('REE - Gt');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,6,6)
% % % %         h=semilogy(C_solid');xlim([1 13]);ylim([1e-3 1e1]);ylabel('C_s^e (ppm)');xlabel('REE - Solid');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %     
% % % %     
% % % %         figure(30);clf(30);
% % % %         subplot(1,5,1)
% % % %         h=semilogy(C_solid');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_s^e (ppm)');xlabel('REE');
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,5,2)
% % % %         h=semilogy(C_total_fluid');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_l^e (ppm)');xlabel('REE'); % instantaneo
% % % %                             hold on;  plot(C_ave','.r', 'MarkerSize',15')
% % % %                             hold on;  plot(C_ave_phi_v','.b', 'MarkerSize',20')
% % % %                             hold on;  plot(C_ave_FM','og', 'MarkerSize',6')
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,5,3)
% % % %         h=semilogy(C_fluid');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_{l,eq}^e (ppm)');xlabel('REE'); % instantaneo
% % % %                             hold on;  plot(C_ave','.r', 'MarkerSize',15')
% % % %                             hold on;  plot(C_ave_phi_v','.b', 'MarkerSize',20')
% % % %                             hold on;  plot(C_ave_FM','og', 'MarkerSize',6')
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,5,4)
% % % %         h=semilogy(C_extracted');xlim([1 13]);ylim([1e-3 1e2]);ylabel('C_{l,ext}^e (ppm)');xlabel('REE'); % instantaneo
% % % %                             hold on;  plot(C_ave','.r', 'MarkerSize',15')
% % % %                             hold on;  plot(C_ave_phi_v','.b', 'MarkerSize',20')
% % % %                             hold on;  plot(C_ave_FM','og', 'MarkerSize',6')
% % % %         set(h, {'color'}, num2cell(Colormap, 2))
% % % %         subplot(1,5,5)
% % % %         plot(100*phi_eq,M_out(:,1)*-1e-3,'.');hold on;plot(100*phi_col,M_out(:,1)*-1e-3,'.');plot(100*F,M_out(:,1)*-1e-3);ylabel('P (kbar)');xlabel('\phi_{eq}, \phi_{ext} and F (%wt)');legend('\phi_{eq}','\phi_{e}','F')
    
    
    %
end

end



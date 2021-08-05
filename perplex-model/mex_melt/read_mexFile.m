%% Adapt the output of the MEX file to the format of the previous calcdata
% this script is adapted to post process more than one minimization. There
% are some known problems 
% a = [a;a];
% b = [b;b];
% c = [c;c];

stringName = cellstr(a);
dataMat    = b;
           
%%%%
% find positions (indices) for relevant variables
log_Ol      = contains(stringName,{'O(JH)'});
log_Cpx     = contains(stringName,{'Cpx(JH)'});
log_Opx     = contains(stringName,{'Opx(JH)'});
log_Gt      = contains(stringName,{'Grt(JH)'});
log_Sp      = contains(stringName,{'Sp(JH)'});
log_Pl      = contains(stringName,{'Pl(JH)'});
log_Melt    = contains(stringName,{'Melt(JH)'});
log_System  = contains(stringName,{'System'});    % allocates space, as "system" is not a variable carried from the mex file
log_Delete  = contains(stringName,{'Delete'});    % allocates space, as "system" is not a variable carried from the mex file
log_System_aux                = find(cellfun('isempty', stringName));
log_aux = log_System_aux(1:end-1) ~= log_System_aux(2:end)-1; log_aux = [true; log_aux];
log_System(log_System_aux(log_aux)) = true;
log_aux = log_System_aux(1:end-1) == log_System_aux(2:end)-1; log_aux = [false; log_aux];
log_Delete(log_System_aux(log_aux)) = true;

% if sum(log_Ol)>1;keyboard;end
% if sum(log_Cpx)>1;keyboard;end
% if sum(log_Opx)>1;keyboard;end
% if sum(log_Gt)>1;keyboard;end
% if sum(log_Sp)>1;keyboard;end
% if sum(log_Pl)>1;keyboard;end
% if sum(log_Melt)>1;keyboard;end


% prepare dataMat in the same format, witgh "System" information at the
% end, and removing the zeros inherited from the mex file.
dataMat(log_System,:) = c;
% dataMat(log_Delete,:) = [];


%%% data format datMat:
% phases:  wt%     vol%  mol%  N(g)  H(J/kg)  S(J/K)  V(J/bar) Cp(J/K)  Alpha(1/K)  Rho(kg/m3)  Vp(km/s)  Vs(km/s)  SiO2(wt%)  Al2O3  FeO  MnO  MgO  CaO  Na2O  Cr2O3  TiO2
% System:  T(K)  P(bar)  G(J)  N(g)  H(J/kg)  S(J/K)  V(J/bar) Cp(J/K)  Alpha(1/K)  Rho(kg/m3)  Vp(km/s)  Vs(km/s)  SiO2(wt%)  Al2O3  FeO  MnO  MgO  CaO  Na2O  Cr2O3  TiO2

%% EXTRACT QUANTITIES OF INTEREST
% entropy, volume, modal composition, densities, thermal expansion, heat capacity

X_Ol = dataMat(log_Ol,1);X_Cpx = dataMat(log_Cpx,1);X_Opx = dataMat(log_Opx,1);
X_Gt = dataMat(log_Gt,1);X_Sp = dataMat(log_Sp,1);X_Pl = dataMat(log_Pl,1);X_Melt = dataMat(log_Melt,1);
X_Ol(isempty(X_Ol)) = 0; X_Cpx(isempty(X_Cpx)) = 0; X_Opx(isempty(X_Opx)) = 0; 
X_Gt(isempty(X_Gt)) = 0; X_Sp(isempty(X_Sp)) = 0; X_Pl(isempty(X_Pl)) = 0; X_Melt(isempty(X_Melt)) = 0;
X_Ol(isnan(X_Ol)) = 0; X_Cpx(isempty(X_Cpx)) = 0; X_Opx(isempty(X_Opx)) = 0; 
X_Gt(isempty(X_Gt)) = 0; X_Sp(isempty(X_Sp)) = 0; X_Pl(isempty(X_Pl)) = 0; X_Melt(isempty(X_Melt)) = 0;
X_Ol(isnan(X_Ol)) = 0; X_Cpx(isnan(X_Cpx)) = 0; X_Opx(isnan(X_Opx)) = 0; 
X_Gt(isnan(X_Gt)) = 0; X_Sp(isnan(X_Sp)) = 0; X_Pl(isnan(X_Pl)) = 0; X_Melt(isnan(X_Melt)) = 0;

Vol_Ol = dataMat(log_Ol,2);Vol_Cpx = dataMat(log_Cpx,2);Vol_Opx = dataMat(log_Opx,2);
Vol_Gt = dataMat(log_Gt,2);Vol_Sp = dataMat(log_Sp,2);Vol_Pl = dataMat(log_Pl,2);Vol_Melt = dataMat(log_Melt,2);
Vol_Ol(isempty(Vol_Ol)) = 0; Vol_Cpx(isempty(Vol_Cpx)) = 0; Vol_Opx(isempty(Vol_Opx)) = 0; 
Vol_Gt(isempty(Vol_Gt)) = 0; Vol_Sp(isempty(Vol_Sp)) = 0; Vol_Pl(isempty(Vol_Pl)) = 0; Vol_Melt(isempty(Vol_Melt)) = 0;
Vol_Ol(isnan(Vol_Ol)) = 0; Vol_Cpx(isnan(Vol_Cpx)) = 0; Vol_Opx(isnan(Vol_Opx)) = 0; 
Vol_Gt(isnan(Vol_Gt)) = 0; Vol_Sp(isnan(Vol_Sp)) = 0; Vol_Pl(isnan(Vol_Pl)) = 0; Vol_Melt(isnan(Vol_Melt)) = 0;

mol_Melt = dataMat(log_Melt,3); mol_Ol = dataMat(log_Ol,3);mol_Cpx = dataMat(log_Cpx,3);mol_Opx = dataMat(log_Opx,3);
mol_Gt = dataMat(log_Gt,3);mol_Sp = dataMat(log_Sp,3);mol_Pl = dataMat(log_Pl,3);
%%% correct for the case that Melt has non-zero wt% but zero mol% (i.e. below precision)
if mol_Melt==0.00; mol_Melt = 0.001; end
mol_Melt(isempty(mol_Melt)) = 0;mol_Ol(isempty(mol_Ol)) = 0;mol_Cpx(isempty(mol_Cpx)) = 0;mol_Opx(isempty(mol_Opx)) = 0;
mol_Gt(isempty(mol_Gt)) = 0;mol_Sp(isempty(mol_Sp)) = 0;mol_Pl(isempty(mol_Pl)) = 0;
mol_Melt(isnan(mol_Melt)) = 0;mol_Ol(isnan(mol_Ol)) = 0;mol_Cpx(isnan(mol_Cpx)) = 0;mol_Opx(isnan(mol_Opx)) = 0;
mol_Gt(isnan(mol_Gt)) = 0;mol_Sp(isnan(mol_Sp)) = 0;mol_Pl(isnan(mol_Pl)) = 0;

N_Melt = dataMat(log_Melt,4); N_Ol = dataMat(log_Ol,4);N_Cpx = dataMat(log_Cpx,4);N_Opx = dataMat(log_Opx,4);
N_Gt = dataMat(log_Gt,4);N_Sp = dataMat(log_Sp,4);N_Pl = dataMat(log_Pl,4);
N_Melt(isempty(N_Melt)) = 0;N_Ol(isempty(N_Ol)) = 0;N_Cpx(isempty(N_Cpx)) = 0;N_Opx(isempty(N_Opx)) = 0;
N_Gt(isempty(N_Gt)) = 0;N_Sp(isempty(N_Sp)) = 0;N_Pl(isempty(N_Pl)) = 0;
N_Melt(isnan(N_Melt)) = 0;N_Ol(isnan(N_Ol)) = 0;N_Cpx(isnan(N_Cpx)) = 0;N_Opx(isnan(N_Opx)) = 0;
N_Gt(isnan(N_Gt)) = 0;N_Sp(isnan(N_Sp)) = 0;N_Pl(isnan(N_Pl)) = 0;

S_Melt = dataMat(log_Melt,6)./(N_Melt*1e-3); S_Ol = dataMat(log_Ol,6)./(N_Ol*1e-3); S_Cpx = dataMat(log_Cpx,6)./(N_Cpx*1e-3); S_Opx = dataMat(log_Opx,6)./(N_Opx*1e-3);
S_Gt = dataMat(log_Gt,6)./(N_Gt*1e-3); S_Sp = dataMat(log_Sp,6)./(N_Sp*1e-3); S_Pl = dataMat(log_Pl,6)./(N_Pl*1e-3);
S_Melt(isempty(S_Melt)) = 0; S_Ol(isempty(S_Ol)) = 0; S_Cpx(isempty(S_Cpx)) = 0; S_Opx(isempty(S_Opx)) = 0; 
S_Gt(isempty(S_Gt)) = 0; S_Sp(isempty(S_Sp)) = 0; S_Pl(isempty(S_Pl)) = 0;
S_Melt(isnan(S_Melt)) = 0; S_Ol(isnan(S_Ol)) = 0; S_Cpx(isnan(S_Cpx)) = 0; S_Opx(isnan(S_Opx)) = 0; 
S_Gt(isnan(S_Gt)) = 0; S_Sp(isnan(S_Sp)) = 0; S_Pl(isnan(S_Pl)) = 0;

V_Melt = dataMat(log_Melt,7)./(N_Melt*1e-3); V_Ol = dataMat(log_Ol,7)./(N_Ol*1e-3); V_Cpx = dataMat(log_Cpx,7)./(N_Cpx*1e-3); V_Opx = dataMat(log_Opx,7)./(N_Opx*1e-3);
V_Gt = dataMat(log_Gt,7)./(N_Gt*1e-3); V_Sp = dataMat(log_Sp,7)./(N_Sp*1e-3); V_Pl = dataMat(log_Pl,7)./(N_Pl*1e-3);
V_Melt(isempty(V_Melt)) = 0; V_Ol(isempty(V_Ol)) = 0; V_Cpx(isempty(V_Cpx)) = 0; V_Opx(isempty(V_Opx)) = 0; 
V_Gt(isempty(V_Gt)) = 0; V_Sp(isempty(V_Sp)) = 0; V_Pl(isempty(V_Pl)) = 0;
V_Melt(isnan(V_Melt)) = 0; V_Ol(isnan(V_Ol)) = 0; V_Cpx(isnan(V_Cpx)) = 0; V_Opx(isnan(V_Opx)) = 0; 
V_Gt(isnan(V_Gt)) = 0; V_Sp(isnan(V_Sp)) = 0; V_Pl(isnan(V_Pl)) = 0;

%%% convert Perplex output Cp [J/K] to "real" Cp [J/(kg K)] by dividing by N [kg];
Cp_Melt = dataMat(log_Melt,8)./(N_Melt*1e-3); Cp_Ol = dataMat(log_Ol,8)./(N_Ol*1e-3);Cp_Cpx = dataMat(log_Cpx,8)./(N_Cpx*1e-3);Cp_Opx = dataMat(log_Opx,8)./(N_Opx*1e-3);
Cp_Gt = dataMat(log_Gt,8)./(N_Gt*1e-3);Cp_Sp = dataMat(log_Sp,8)./(N_Sp*1e-3);Cp_Pl = dataMat(log_Pl,8)./(N_Pl*1e-3);
Cp_Melt(isempty(Cp_Melt)) = 0;Cp_Ol(isempty(Cp_Ol)) = 0; Cp_Cpx(isempty(Cp_Cpx)) = 0; Cp_Opx(isempty(Cp_Opx)) = 0; 
Cp_Gt(isempty(Cp_Gt)) = 0; Cp_Sp(isempty(Cp_Sp)) = 0; Cp_Pl(isempty(Cp_Pl)) = 0;
Cp_Melt(isnan(Cp_Melt)) = 0;Cp_Ol(isnan(Cp_Ol)) = 0; Cp_Cpx(isnan(Cp_Cpx)) = 0; Cp_Opx(isnan(Cp_Opx)) = 0; 
Cp_Gt(isnan(Cp_Gt)) = 0; Cp_Sp(isnan(Cp_Sp)) = 0; Cp_Pl(isnan(Cp_Pl)) = 0;

% alpha_Melt = dataMat(log_Melt,9)./(N_Melt*1e-3); alpha_Ol = dataMat(log_Ol,9)./(N_Ol*1e-3);alpha_Cpx = dataMat(log_Cpx,9)./(N_Cpx*1e-3);alpha_Opx = dataMat(log_Opx,9)./(N_Opx*1e-3);
% alpha_Gt = dataMat(log_Gt,9)./(N_Gt*1e-3);alpha_Sp = dataMat(log_Sp,9)./(N_Sp*1e-3);alpha_Pl = dataMat(log_Pl,9)./(N_Pl*1e-3);
% alpha_Melt(isempty(alpha_Melt)) = 0; alpha_Ol(isempty(alpha_Ol)) = 0; alpha_Cpx(isempty(alpha_Cpx)) = 0; alpha_Opx(isempty(alpha_Opx)) = 0; 
% alpha_Gt(isempty(alpha_Gt)) = 0; alpha_Sp(isempty(alpha_Sp)) = 0; alpha_Pl(isempty(alpha_Pl)) = 0;

alpha_Melt = dataMat(log_Melt,9); alpha_Ol = dataMat(log_Ol,9);alpha_Cpx = dataMat(log_Cpx,9);alpha_Opx = dataMat(log_Opx,9);
alpha_Gt = dataMat(log_Gt,9);alpha_Sp = dataMat(log_Sp,9);alpha_Pl = dataMat(log_Pl,9);
alpha_Melt(isempty(alpha_Melt)) = 0; alpha_Ol(isempty(alpha_Ol)) = 0; alpha_Cpx(isempty(alpha_Cpx)) = 0; alpha_Opx(isempty(alpha_Opx)) = 0; 
alpha_Gt(isempty(alpha_Gt)) = 0; alpha_Sp(isempty(alpha_Sp)) = 0; alpha_Pl(isempty(alpha_Pl)) = 0;
alpha_Melt(isnan(alpha_Melt)) = 0; alpha_Ol(isnan(alpha_Ol)) = 0; alpha_Cpx(isnan(alpha_Cpx)) = 0; alpha_Opx(isnan(alpha_Opx)) = 0; 
alpha_Gt(isnan(alpha_Gt)) = 0; alpha_Sp(isnan(alpha_Sp)) = 0; alpha_Pl(isnan(alpha_Pl)) = 0;

Rho_Ol = dataMat(log_Ol,10);Rho_Cpx = dataMat(log_Cpx,10);Rho_Opx = dataMat(log_Opx,10);
Rho_Gt = dataMat(log_Gt,10);Rho_Sp = dataMat(log_Sp,10);Rho_Pl = dataMat(log_Pl,10);Rho_Melt = dataMat(log_Melt,10);
Rho_Ol(isempty(Rho_Ol)) = 0; Rho_Cpx(isempty(Rho_Cpx)) = 0;Rho_Opx(isempty(Rho_Opx)) = 0;
Rho_Gt(isempty(Rho_Gt)) = 0;Rho_Sp(isempty(Rho_Sp)) = 0;Rho_Pl(isempty(Rho_Pl)) = 0;Rho_Melt(isempty(Rho_Melt)) = 0;
Rho_Ol(isnan(Rho_Ol)) = 0; Rho_Cpx(isnan(Rho_Cpx)) = 0;Rho_Opx(isnan(Rho_Opx)) = 0;
Rho_Gt(isnan(Rho_Gt)) = 0;Rho_Sp(isnan(Rho_Sp)) = 0;Rho_Pl(isnan(Rho_Pl)) = 0;Rho_Melt(isnan(Rho_Melt)) = 0;

Vp_Ol = dataMat(log_Ol,11);Vp_Cpx = dataMat(log_Cpx,11);Vp_Opx = dataMat(log_Opx,11);
Vp_Gt = dataMat(log_Gt,11);Vp_Sp = dataMat(log_Sp,11);Vp_Pl = dataMat(log_Pl,11);Vp_Melt = dataMat(log_Melt,11);
Vp_Ol(isempty(Vp_Ol)) = 0; Vp_Cpx(isempty(Vp_Cpx)) = 0;Vp_Opx(isempty(Vp_Opx)) = 0;
Vp_Gt(isempty(Vp_Gt)) = 0;Vp_Sp(isempty(Vp_Sp)) = 0;Vp_Pl(isempty(Vp_Pl)) = 0;Vp_Melt(isempty(Vp_Melt)) = 0;
Vp_Ol(isnan(Vp_Ol)) = 0; Vp_Cpx(isnan(Vp_Cpx)) = 0;Vp_Opx(isnan(Vp_Opx)) = 0;
Vp_Gt(isnan(Vp_Gt)) = 0;Vp_Sp(isnan(Vp_Sp)) = 0;Vp_Pl(isnan(Vp_Pl)) = 0;Vp_Melt(isnan(Vp_Melt)) = 0;

Vs_Ol = dataMat(log_Ol,12);Vs_Cpx = dataMat(log_Cpx,12);Vs_Opx = dataMat(log_Opx,12);
Vs_Gt = dataMat(log_Gt,12);Vs_Sp = dataMat(log_Sp,12);Vs_Pl = dataMat(log_Pl,12);Vs_Melt = dataMat(log_Melt,12);
Vs_Ol(isempty(Vs_Ol)) = 0; Vs_Cpx(isempty(Vs_Cpx)) = 0;Vs_Opx(isempty(Vs_Opx)) = 0;
Vs_Gt(isempty(Vs_Gt)) = 0;Vs_Sp(isempty(Vs_Sp)) = 0;Vs_Pl(isempty(Vs_Pl)) = 0;Vs_Melt(isempty(Vs_Melt)) = 0;
Vs_Ol(isnan(Vs_Ol)) = 0; Vs_Cpx(isnan(Vs_Cpx)) = 0;Vs_Opx(isnan(Vs_Opx)) = 0;
Vs_Gt(isnan(Vs_Gt)) = 0;Vs_Sp(isnan(Vs_Sp)) = 0;Vs_Pl(isnan(Vs_Pl)) = 0;Vs_Melt(isnan(Vs_Melt)) = 0;

H_Ol = dataMat(log_Ol,5)./(N_Ol*1e-3); H_Cpx = dataMat(log_Cpx,5)./(N_Cpx*1e-3); H_Opx = dataMat(log_Opx,5)./(N_Opx*1e-3);
H_Gt = dataMat(log_Gt,5)./(N_Gt*1e-3); H_Sp = dataMat(log_Sp,5)./(N_Sp*1e-3); H_Pl = dataMat(log_Pl,5)./(N_Pl*1e-3); H_Melt = dataMat(log_Melt,5)./(N_Melt*1e-3);
H_Ol(isempty(H_Ol)) = 0; H_Cpx(isempty(H_Cpx)) = 0; H_Opx(isempty(H_Opx)) = 0;
H_Gt(isempty(H_Gt)) = 0;H_Sp(isempty(H_Sp)) = 0;H_Pl(isempty(H_Pl)) = 0;H_Melt(isempty(H_Melt)) = 0;
H_Ol(isnan(H_Ol)) = 0; H_Cpx(isnan(H_Cpx)) = 0; H_Opx(isnan(H_Opx)) = 0;
H_Gt(isnan(H_Gt)) = 0;H_Sp(isnan(H_Sp)) = 0;H_Pl(isnan(H_Pl)) = 0;H_Melt(isnan(H_Melt)) = 0;

C_Melt = dataMat(log_Melt,13:21);C_Ol = dataMat(log_Ol,13:21);C_Cpx = dataMat(log_Cpx,13:21);C_Opx = dataMat(log_Opx,13:21);
C_Gt = dataMat(log_Gt,13:21);C_Sp = dataMat(log_Sp,13:21);C_Pl = dataMat(log_Pl,13:21);
C_Melt(isempty(C_Melt)) = 0;C_Ol(isempty(C_Ol)) = 0;C_Cpx(isempty(C_Cpx)) = 0;C_Opx(isempty(C_Opx)) = 0;
C_Gt(isempty(C_Gt)) = 0;C_Sp(isempty(C_Sp)) = 0;C_Pl(isempty(C_Pl)) = 0;
C_Melt(isnan(C_Melt)) = 0;C_Ol(isnan(C_Ol)) = 0;C_Cpx(isnan(C_Cpx)) = 0;C_Opx(isnan(C_Opx)) = 0;
C_Gt(isnan(C_Gt)) = 0;C_Sp(isnan(C_Sp)) = 0;C_Pl(isnan(C_Pl)) = 0;

%%% correct for duplicate phases in Perple_X output
if length(X_Sp) > 1 
Rho_Sp = sum((Rho_Sp.*X_Sp))/sum(X_Sp); Vp_Sp = sum((Vp_Sp.*X_Sp))/sum(X_Sp); Vs_Sp = sum((Vs_Sp.*X_Sp))/sum(X_Sp); H_Sp = sum((H_Sp.*X_Sp))/sum(X_Sp); N_Sp = sum(N_Sp); C_Sp = sum((C_Sp.*X_Sp))/sum(X_Sp);
V_Sp = sum((V_Sp.*X_Sp))/sum(X_Sp); S_Sp = sum((S_Sp.*X_Sp))/sum(X_Sp); alpha_Sp = sum((alpha_Sp.*X_Sp))/sum(X_Sp); Cp_Sp = sum((Cp_Sp.*mol_Sp))/sum(mol_Sp);
X_Sp = sum(X_Sp); Vol_Sp = sum(Vol_Sp); mol_Sp = sum(mol_Sp); 
end
Rho_Sp(isnan(Rho_Sp)) = 0; Vp_Sp(isnan(Vp_Sp)) = 0; Vs_Sp(isnan(Vs_Sp)) = 0; H_Sp(isnan(H_Sp)) = 0; N_Sp(isnan(N_Sp)) = 0; C_Sp(isnan(C_Sp)) = 0; 
V_Sp(isnan(V_Sp)) = 0; S_Sp(isnan(S_Sp)) = 0; alpha_Sp(isnan(alpha_Sp)) = 0; Cp_Sp(isnan(Cp_Sp)) = 0; 
X_Sp(isnan(X_Sp)) = 0; Vol_Sp(isnan(Vol_Sp)) = 0; mol_Sp(isnan(mol_Sp)) = 0; 

if length(X_Ol) > 1 
Rho_Ol = sum((Rho_Ol.*X_Ol))/sum(X_Ol); Vp_Ol = sum((Vp_Ol.*X_Ol))/sum(X_Ol); Vs_Ol = sum((Vs_Ol.*X_Ol))/sum(X_Ol); H_Ol = sum((H_Ol.*X_Ol))/sum(X_Ol); N_Ol = sum(N_Ol); C_Ol = sum((C_Ol.*X_Ol))/sum(X_Ol);
V_Ol = sum((V_Ol.*X_Ol))/sum(X_Ol); S_Ol = sum((S_Ol.*X_Ol))/sum(X_Ol); alpha_Ol = sum((alpha_Ol.*X_Ol))/sum(X_Ol); Cp_Ol = sum((Cp_Ol.*mol_Ol))/sum(mol_Ol);
X_Ol = sum(X_Ol);Vol_Ol = sum(Vol_Ol); mol_Ol = sum(mol_Ol); 
end 
Rho_Ol(isnan(Rho_Ol)) = 0; Vp_Ol(isnan(Vp_Ol)) = 0; Vs_Ol(isnan(Vs_Ol)) = 0; H_Ol(isnan(H_Ol)) = 0; N_Ol(isnan(N_Ol)) = 0; C_Ol(isnan(C_Ol)) = 0; 
V_Ol(isnan(V_Ol)) = 0; S_Ol(isnan(S_Ol)) = 0; alpha_Ol(isnan(alpha_Ol)) = 0; Cp_Ol(isnan(Cp_Ol)) = 0; 
X_Ol(isnan(X_Ol)) = 0; Vol_Ol(isnan(Vol_Ol)) = 0; mol_Ol(isnan(mol_Ol)) = 0; 

if length(X_Cpx) > 1 
Rho_Cpx = sum((Rho_Cpx.*X_Cpx))/sum(X_Cpx); Vp_Cpx = sum((Vp_Cpx.*X_Cpx))/sum(X_Cpx); Vs_Cpx = sum((Vs_Cpx.*X_Cpx))/sum(X_Cpx); H_Cpx = sum((H_Cpx.*X_Cpx))/sum(X_Cpx); N_Cpx = sum(N_Cpx); C_Cpx = sum((C_Cpx.*X_Cpx))/sum(X_Cpx);
V_Cpx = sum((V_Cpx.*X_Cpx))/sum(X_Cpx); S_Cpx = sum((S_Cpx.*X_Cpx))/sum(X_Cpx); alpha_Cpx = sum((alpha_Cpx.*mol_Cpx))/sum(mol_Cpx); Cp_Cpx = sum((Cp_Cpx.*X_Cpx))/sum(X_Cpx);
X_Cpx = sum(X_Cpx); Vol_Cpx = sum(Vol_Cpx); mol_Cpx = sum(mol_Cpx);  
end
Rho_Cpx(isnan(Rho_Cpx)) = 0; Vp_Cpx(isnan(Vp_Cpx)) = 0; Vs_Cpx(isnan(Vs_Cpx)) = 0; H_Cpx(isnan(H_Cpx)) = 0; N_Cpx(isnan(N_Cpx)) = 0; C_Cpx(isnan(C_Cpx)) = 0; 
V_Cpx(isnan(V_Cpx)) = 0; S_Cpx(isnan(S_Cpx)) = 0; alpha_Cpx(isnan(alpha_Cpx)) = 0; Cp_Cpx(isnan(Cp_Cpx)) = 0; 
X_Cpx(isnan(X_Cpx)) = 0; Vol_Cpx(isnan(Vol_Cpx)) = 0; mol_Cpx(isnan(mol_Cpx)) = 0; 

if length(X_Opx) > 1 
Rho_Opx = sum((Rho_Opx.*X_Opx))/sum(X_Opx); Vp_Opx = sum((Vp_Opx.*X_Opx))/sum(X_Opx); Vs_Opx = sum((Vs_Opx.*X_Opx))/sum(X_Opx); H_Opx = sum((H_Opx.*X_Opx))/sum(X_Opx); N_Opx = sum(N_Opx); C_Opx = sum((C_Opx.*X_Opx))/sum(X_Opx);
V_Opx = sum((V_Opx.*X_Opx))/sum(X_Opx); S_Opx = sum((S_Opx.*X_Opx))/sum(X_Opx); alpha_Opx = sum((alpha_Opx.*X_Opx))/sum(X_Opx); Cp_Opx = sum((Cp_Opx.*mol_Opx))/sum(mol_Opx);
X_Opx = sum(X_Opx); Vol_Opx = sum(Vol_Opx); mol_Opx = sum(mol_Opx); 
end 
Rho_Opx(isnan(Rho_Opx)) = 0; Vp_Opx(isnan(Vp_Opx)) = 0; Vs_Opx(isnan(Vs_Opx)) = 0; H_Opx(isnan(H_Opx)) = 0; N_Opx(isnan(N_Opx)) = 0; C_Opx(isnan(C_Opx)) = 0; 
V_Opx(isnan(V_Opx)) = 0; S_Opx(isnan(S_Opx)) = 0; alpha_Opx(isnan(alpha_Opx)) = 0; Cp_Opx(isnan(Cp_Opx)) = 0; 
X_Opx(isnan(X_Opx)) = 0; Vol_Opx(isnan(Vol_Opx)) = 0; mol_Opx(isnan(mol_Opx)) = 0; 

if length(X_Gt) > 1 
Rho_Gt = sum((Rho_Gt.*X_Gt))/sum(X_Gt); Vp_Gt = sum((Vp_Gt.*X_Gt))/sum(X_Gt); Vs_Gt = sum((Vs_Gt.*X_Gt))/sum(X_Gt); H_Gt = sum((H_Gt.*X_Gt))/sum(X_Gt); N_Gt = sum(N_Gt); C_Gt = sum((C_Gt.*X_Gt))/sum(X_Gt);
V_Gt = sum((V_Gt.*X_Gt))/sum(X_Gt); S_Gt = sum((S_Gt.*X_Gt))/sum(X_Gt); alpha_Gt = sum((alpha_Gt.*X_Gt))/sum(X_Gt); Cp_Gt = sum((Cp_Gt.*mol_Gt))/sum(mol_Gt);
X_Gt = sum(X_Gt); Vol_Gt = sum(Vol_Gt); mol_Gt = sum(mol_Gt);  
end
Rho_Gt(isnan(Rho_Gt)) = 0; Vp_Gt(isnan(Vp_Gt)) = 0; Vs_Gt(isnan(Vs_Gt)) = 0; H_Gt(isnan(H_Gt)) = 0; N_Gt(isnan(N_Gt)) = 0; C_Gt(isnan(C_Gt)) = 0; 
V_Gt(isnan(V_Gt)) = 0; S_Gt(isnan(S_Gt)) = 0; alpha_Gt(isnan(alpha_Gt)) = 0; Cp_Gt(isnan(Cp_Gt)) = 0; 
X_Gt(isnan(X_Gt)) = 0; Vol_Gt(isnan(Vol_Gt)) = 0; mol_Gt(isnan(mol_Gt)) = 0; 

if length(X_Pl) > 1 
Rho_Pl = sum((Rho_Pl.*X_Pl))/sum(X_Pl); Vp_Pl = sum((Vp_Pl.*X_Pl))/sum(X_Pl); Vs_Pl = sum((Vs_Pl.*X_Pl))/sum(X_Pl); H_Pl = sum((H_Pl.*X_Pl))/sum(X_Pl); N_Pl = sum(N_Pl); C_Pl = sum((C_Pl.*X_Pl))/sum(X_Pl);
V_Pl = sum((V_Pl.*X_Pl))/sum(X_Pl); S_Pl = sum((S_Pl.*X_Pl))/sum(X_Pl); alpha_Pl = sum((alpha_Pl.*X_Pl))/sum(X_Pl); Cp_Pl = sum((Cp_Pl.*mol_Pl))/sum(mol_Pl);
X_Pl = sum(X_Pl); Vol_Pl = sum(Vol_Pl); mol_Pl = sum(mol_Pl);  
end
Rho_Pl(isnan(Rho_Pl)) = 0; Vp_Pl(isnan(Vp_Pl)) = 0; Vs_Pl(isnan(Vs_Pl)) = 0; H_Pl(isnan(H_Pl)) = 0; N_Pl(isnan(N_Pl)) = 0; C_Pl(isnan(C_Pl)) = 0; 
V_Pl(isnan(V_Pl)) = 0; S_Pl(isnan(S_Pl)) = 0; alpha_Pl(isnan(alpha_Pl)) = 0; Cp_Pl(isnan(Cp_Pl)) = 0; 
X_Pl(isnan(X_Pl)) = 0; Vol_Pl(isnan(Vol_Pl)) = 0; mol_Pl(isnan(mol_Pl)) = 0; 

if length(X_Melt) > 1 
Rho_Melt = sum((Rho_Melt.*X_Melt))/sum(X_Melt); Vp_Melt = sum((Vp_Melt.*X_Melt))/sum(X_Melt); Vs_Melt = sum((Vs_Melt.*X_Melt))/sum(X_Melt); H_Melt = sum((H_Melt.*X_Melt))/sum(X_Melt); N_Melt = sum(N_Melt); C_Melt = sum((C_Melt.*X_Melt))/sum(X_Melt);
V_Melt = sum((V_Melt.*X_Melt))/sum(X_Melt); S_Melt = sum((S_Melt.*X_Melt))/sum(X_Melt); alpha_Melt = sum((alpha_Melt.*X_Melt))/sum(X_Melt); Cp_Melt = sum((Cp_Melt.*mol_Melt))/sum(mol_Melt);
X_Melt = sum(X_Melt); Vol_Melt = sum(Vol_Melt); mol_Melt = sum(mol_Melt);  
end
Rho_Melt(isnan(Rho_Melt)) = 0; Vp_Melt(isnan(Vp_Melt)) = 0; Vs_Melt(isnan(Vs_Melt)) = 0; H_Melt(isnan(H_Melt)) = 0; N_Melt(isnan(N_Melt)) = 0; C_Melt(isnan(C_Melt)) = 0; 
V_Melt(isnan(V_Melt)) = 0; S_Melt(isnan(S_Melt)) = 0; alpha_Melt(isnan(alpha_Melt)) = 0; Cp_Melt(isnan(Cp_Melt)) = 0; 
X_Melt(isnan(X_Melt)) = 0; Vol_Melt(isnan(Vol_Melt)) = 0; mol_Melt(isnan(mol_Melt)) = 0; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate System values as (normalised) sum of 7 phases

%N_tot = N_Melt + N_Ol+ N_Cpx+ N_Opx+ N_Gt+ N_Sp+ N_Pl;

mol_tot = mol_Melt + mol_Ol+ mol_Cpx+ mol_Opx+ mol_Gt+ mol_Sp+ mol_Pl;
Mol_Melt = mol_Melt/mol_tot*100; Mol_Ol = mol_Ol/mol_tot*100; Mol_Cpx = mol_Cpx/mol_tot*100; Mol_Opx = mol_Opx/mol_tot*100; 
Mol_Gt = mol_Gt/mol_tot*100; Mol_Sp = mol_Sp/mol_tot*100; Mol_Pl = mol_Pl/mol_tot*100; 
Mol_System = Mol_Melt + Mol_Ol+ Mol_Cpx+ Mol_Opx+ Mol_Gt+ Mol_Sp+ Mol_Pl;

X_tot = X_Melt + X_Ol+ X_Cpx+ X_Opx+ X_Gt+ X_Sp+ X_Pl;
Xn_Melt = X_Melt/X_tot*100; Xn_Ol = X_Ol/X_tot*100; Xn_Cpx = X_Cpx/X_tot*100; Xn_Opx = X_Opx/X_tot*100; 
Xn_Gt = X_Gt/X_tot*100; Xn_Sp = X_Sp/X_tot*100; Xn_Pl = X_Pl/X_tot*100; 
Xn_System = Xn_Melt + Xn_Ol+ Xn_Cpx+ Xn_Opx+ Xn_Gt+ Xn_Sp+ Xn_Pl;

Vol_tot = Vol_Melt + Vol_Ol+ Vol_Cpx+ Vol_Opx+ Vol_Gt+ Vol_Sp+ Vol_Pl;
Voln_Melt = Vol_Melt/Vol_tot*100; Voln_Ol = Vol_Ol/Vol_tot*100; Voln_Cpx = Vol_Cpx/Vol_tot*100; Voln_Opx = Vol_Opx/Vol_tot*100; 
Voln_Gt = Vol_Gt/Vol_tot*100; Voln_Sp = Vol_Sp/Vol_tot*100; Voln_Pl = Vol_Pl/Vol_tot*100; 
Voln_System = Voln_Melt + Voln_Ol+ Voln_Cpx+ Voln_Opx+ Voln_Gt+ Voln_Sp+ Voln_Pl;

Cn_Melt = C_Melt*Xn_Melt/Xn_System; Cn_Ol = C_Ol*Xn_Ol/Xn_System; Cn_Cpx = C_Cpx*Xn_Cpx/Xn_System; Cn_Opx = C_Opx*Xn_Opx/Xn_System;
Cn_Gt = C_Gt*Xn_Gt/Xn_System; Cn_Sp = C_Sp*Xn_Sp/Xn_System; Cn_Pl = C_Pl*Xn_Pl/Xn_System; 

C_System =  (Cn_Melt + Cn_Ol+ Cn_Cpx+ Cn_Opx+ Cn_Gt+ Cn_Sp+ Cn_Pl);

% if round(sum(C_System))~=100
%     error('System not normalised to 100')
% end

P = dataMat(log_System,2); T = dataMat(log_System,1);

Cp_System   = (Cp_Melt*Mol_Melt + Cp_Ol*Mol_Ol+ Cp_Cpx*Mol_Cpx+ Cp_Opx*Mol_Opx+ Cp_Gt*Mol_Gt+ Cp_Sp*Mol_Sp+ Cp_Pl*Mol_Pl) / Mol_System;

S_System    = (S_Melt*Xn_Melt + S_Ol*Xn_Ol+ S_Cpx*Xn_Cpx+ S_Opx*Xn_Opx+ S_Gt*Xn_Gt+ S_Sp*Xn_Sp+ S_Pl*Xn_Pl) / Xn_System;
alpha_System = (alpha_Melt*Xn_Melt + alpha_Ol*Xn_Ol+ alpha_Cpx*Xn_Cpx+ alpha_Opx*Xn_Opx+ alpha_Gt*Xn_Gt+ alpha_Sp*Xn_Sp+ alpha_Pl*Xn_Pl) / Xn_System;
V_System    = (V_Melt*Xn_Melt + V_Ol*Xn_Ol+ V_Cpx*Xn_Cpx+ V_Opx*Xn_Opx+ V_Gt*Xn_Gt+ V_Sp*Xn_Sp+ V_Pl*Xn_Pl) / Xn_System;
rho_System  = (Rho_Melt*Xn_Melt + Rho_Ol*Xn_Ol+ Rho_Cpx*Xn_Cpx+ Rho_Opx*Xn_Opx+ Rho_Gt*Xn_Gt+ Rho_Sp*Xn_Sp+ Rho_Pl*Xn_Pl) / Xn_System;
Vp_System  = (Vp_Melt*Xn_Melt + Vp_Ol*Xn_Ol+ Vp_Cpx*Xn_Cpx+ Vp_Opx*Xn_Opx+ Vp_Gt*Xn_Gt+ Vp_Sp*Xn_Sp+ Vp_Pl*Xn_Pl) / Xn_System;
Vs_System  = (Vs_Melt*Xn_Melt + Vs_Ol*Xn_Ol+ Vs_Cpx*Xn_Cpx+ Vs_Opx*Xn_Opx+ Vs_Gt*Xn_Gt+ Vs_Sp*Xn_Sp+ Vs_Pl*Xn_Pl) / Xn_System;
H_System    = (H_Melt*Xn_Melt + H_Ol*Xn_Ol+ H_Cpx*Xn_Cpx+ H_Opx*Xn_Opx+ H_Gt*Xn_Gt+ H_Sp*Xn_Sp+ H_Pl*Xn_Pl) / Xn_System;
H_System(isempty(H_System)) = 0;

% Solid residue
% phases:  wt%     vol%  mol%  N(g)  H(J/kg)  S(J/K)  V(J/bar) Cp(J/K)  Alpha(1/K)  Rho(kg/m3)  SiO2(wt%)  Al2O3  FeO  MnO  MgO  CaO  Na2O  Cr2O3  TiO2  [K2O]
% System:  T(K)  P(bar)  G(J)  N(g)  H(J/kg)  S(J/K)  V(J/bar) Cp(J/K)  Alpha(1/K)  Rho(kg/m3)  SiO2(wt%)  Al2O3  FeO  MnO  MgO  CaO  Na2O  Cr2O3  TiO2  [K2O]

C_Solid = (X_Ol*C_Ol + X_Cpx*C_Cpx + X_Opx*C_Opx + X_Gt*C_Gt + X_Sp*C_Sp + X_Pl*C_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);
X_Solid = (X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl)/(X_Melt + X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);
Vol_Solid = (Vol_Ol + Vol_Cpx + Vol_Opx + Vol_Gt + Vol_Sp + Vol_Pl)/(Vol_Melt + Vol_Ol + Vol_Cpx + Vol_Opx + Vol_Gt + Vol_Sp + Vol_Pl);
Mol_Solid = (mol_Ol + mol_Cpx + mol_Opx + mol_Gt + mol_Sp + mol_Pl)/mol_tot*100;
H_Solid    = (H_Ol*X_Ol+ H_Cpx*X_Cpx+ H_Opx*X_Opx+ H_Gt*X_Gt+ H_Sp*X_Sp+ H_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);    % J/kg
S_Solid    = (S_Ol*X_Ol+ S_Cpx*X_Cpx+ S_Opx*X_Opx+ S_Gt*X_Gt+ S_Sp*X_Sp+ S_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);    % J/kg/K 
V_Solid    = (V_Ol*X_Ol+ V_Cpx*X_Cpx+ V_Opx*X_Opx+ V_Gt*X_Gt+ V_Sp*X_Sp+ V_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);
% Cp_Solid   = (Cp_Ol*Mol_Ol+ Cp_Cpx*Mol_Cpx+ Cp_Opx*Mol_Opx+ Cp_Gt*Mol_Gt+ Cp_Sp*Mol_Sp+ Cp_Pl*Mol_Pl)/(mol_Ol + mol_Cpx + mol_Opx + mol_Gt + mol_Sp + mol_Pl);
Cp_Solid   = (Cp_Ol*X_Ol+ Cp_Cpx*X_Cpx+ Cp_Opx*X_Opx+ Cp_Gt*X_Gt+ Cp_Sp*X_Sp+ Cp_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);  % J/kg/K 
alpha_Solid = (alpha_Ol*X_Ol+ alpha_Cpx*X_Cpx+ alpha_Opx*X_Opx+ alpha_Gt*X_Gt+ alpha_Sp*X_Sp+ alpha_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);
rho_Solid  = (Rho_Ol*Vol_Ol+ Rho_Cpx*Vol_Cpx+ Rho_Opx*Vol_Opx+ Rho_Gt*Vol_Gt+ Rho_Sp*Vol_Sp+ Rho_Pl*Vol_Pl)/(Vol_Ol + Vol_Cpx + Vol_Opx + Vol_Gt + Vol_Sp + Vol_Pl);
Vp_Solid = (Vp_Ol*X_Ol+ Vp_Cpx*X_Cpx+ Vp_Opx*X_Opx+ Vp_Gt*X_Gt+ Vp_Sp*X_Sp+ Vp_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);
Vs_Solid = (Vs_Ol*X_Ol+ Vs_Cpx*X_Cpx+ Vs_Opx*X_Opx+ Vs_Gt*X_Gt+ Vs_Sp*X_Sp+ Vs_Pl*X_Pl)/(X_Ol + X_Cpx + X_Opx + X_Gt + X_Sp + X_Pl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% OLD
% 
% % entropy, volume, modal composition, densities, thermal expansion, heat capacity
% P = dataMat(logical_System,1); T = dataMat(logical_System,2);
% S = dataMat(logical_System,6); Cp_System = dataMat(logical_System,8); alpha_System = dataMat(logical_System,9); V_System = dataMat(logical_System,7); rho_System = dataMat(logical_System,10);
% 
% S_Melt = dataMat(logical_Melt,6); S_Melt(isempty(S_Melt)) = 0;
% V_Melt = dataMat(logical_Melt,7); V_Melt(isempty(V_Melt)) = 0;
% 
% X_Ol = dataMat(logical_Ol,1);X_Cpx = dataMat(logical_Cpx,1);X_Opx = dataMat(logical_Opx,1);X_Gt = dataMat(logical_Gt,1);X_Sp = dataMat(logical_Sp,1);X_Pl = dataMat(logical_Pl,1);X_Melt = dataMat(logical_Melt,1);
% X_Ol(isempty(X_Ol)) = 0; X_Cpx(isempty(X_Cpx)) = 0; X_Opx(isempty(X_Opx)) = 0; X_Gt(isempty(X_Gt)) = 0; X_Sp(isempty(X_Sp)) = 0; X_Pl(isempty(X_Pl)) = 0; X_Melt(isempty(X_Melt)) = 0;
% 
% Rho_Ol = dataMat(logical_Ol,10);Rho_Cpx = dataMat(logical_Cpx,10);Rho_Opx = dataMat(logical_Opx,10);Rho_Gt = dataMat(logical_Gt,10);Rho_Sp = dataMat(logical_Sp,10);Rho_Pl = dataMat(logical_Pl,10);Rho_Melt = dataMat(logical_Melt,10);
% Rho_Ol(isempty(Rho_Ol)) = 0; Rho_Cpx(isempty(Rho_Cpx)) = 0;Rho_Opx(isempty(Rho_Opx)) = 0;Rho_Gt(isempty(Rho_Gt)) = 0;Rho_Sp(isempty(Rho_Sp)) = 0;Rho_Pl(isempty(Rho_Pl)) = 0;Rho_Melt(isempty(Rho_Melt)) = 0;
% 
% N_Melt = dataMat(logical_Melt,4); N_Melt(isempty(N_Melt)) = 0; N_System = dataMat(logical_System,4);
% 
% H_Ol = dataMat(logical_Ol,5);H_Cpx = dataMat(logical_Cpx,5);H_Opx = dataMat(logical_Opx,5);H_Gt = dataMat(logical_Gt,5);H_Sp = dataMat(logical_Sp,5);H_Pl = dataMat(logical_Pl,5);H_Melt = dataMat(logical_Melt,5);H_System = dataMat(logical_System,5);
% H_Ol(isempty(H_Ol)) = 0; H_Cpx(isempty(H_Cpx)) = 0;H_Opx(isempty(H_Opx)) = 0;H_Gt(isempty(H_Gt)) = 0;H_Sp(isempty(H_Sp)) = 0;H_Pl(isempty(H_Pl)) = 0;H_Melt(isempty(H_Melt)) = 0;H_System(isempty(H_System)) = 0;
% 
% %%% correct for double phases in Perple_X output
% X_Sp = mean(X_Sp); Rho_Sp = mean(Rho_Sp); H_Sp = mean(H_Sp);


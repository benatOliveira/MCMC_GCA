function [P, T, S_System, S_Solid, S_Melt, V_Solid, V_Melt, Cp_Solid, alpha_Solid, rho_Solid,...
    X_Melt, X_Ol, X_Cpx, X_Opx, X_Gt, X_Sp, X_Pl, Vol_Melt, Vol_Ol, Vol_Cpx, Vol_Opx, Vol_Gt, Vol_Sp, Vol_Pl,...
    Rho_Melt, Rho_Ol, Rho_Cpx, Rho_Opx, Rho_Gt, Rho_Sp, Rho_Pl,...
    H_Solid, H_Melt, H_Gt, H_Sp,... % 1/dFdT_P dTdP_F];
    C_Solid, C_Melt, C_Cpx, C_Gt] = ...
    interpolatePerplex(P_par, T_par, S_System_par, S_Solid_par, S_Melt_par, V_Solid_par, V_Melt_par, Cp_Solid_par, alpha_Solid_par, rho_Solid_par,...
    X_Melt_par, X_Ol_par, X_Cpx_par, X_Opx_par, X_Gt_par, X_Sp_par, X_Pl_par, Vol_Melt_par, Vol_Ol_par, Vol_Cpx_par, Vol_Opx_par, Vol_Gt_par, Vol_Sp_par, Vol_Pl_par,...
    Rho_Melt_par, Rho_Ol_par, Rho_Cpx_par, Rho_Opx_par, Rho_Gt_par, Rho_Sp_par, Rho_Pl_par,...
    H_Solid_par, H_Melt_par, H_Gt_par, H_Sp_par,... % 1/dFdT_P dTdP_F];
    C_Solid_par, C_Melt_par, C_Cpx_par, C_Gt_par,DeltaNa2O,index_Na2O)


P = P_par{index_Na2O(1)}*DeltaNa2O + P_par{index_Na2O(2)}*(1-DeltaNa2O);
T = T_par{index_Na2O(1)}*DeltaNa2O + T_par{index_Na2O(2)}*(1-DeltaNa2O);

S_System = S_System_par{index_Na2O(1)}*DeltaNa2O + S_System_par{index_Na2O(2)}*(1-DeltaNa2O);
S_Solid = S_Solid_par{index_Na2O(1)}*DeltaNa2O + S_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);
S_Melt = S_Melt_par{index_Na2O(1)}*DeltaNa2O + S_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
V_Solid = V_Solid_par{index_Na2O(1)}*DeltaNa2O + V_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);
V_Melt = V_Melt_par{index_Na2O(1)}*DeltaNa2O + V_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
Cp_Solid = Cp_Solid_par{index_Na2O(1)}*DeltaNa2O + Cp_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);
alpha_Solid = alpha_Solid_par{index_Na2O(1)}*DeltaNa2O + alpha_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);
rho_Solid = rho_Solid_par{index_Na2O(1)}*DeltaNa2O + rho_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);


X_Melt = X_Melt_par{index_Na2O(1)}*DeltaNa2O + X_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
X_Ol = X_Ol_par{index_Na2O(1)}*DeltaNa2O + X_Ol_par{index_Na2O(2)}*(1-DeltaNa2O);
X_Cpx = X_Cpx_par{index_Na2O(1)}*DeltaNa2O + X_Cpx_par{index_Na2O(2)}*(1-DeltaNa2O);
X_Opx = X_Opx_par{index_Na2O(1)}*DeltaNa2O + X_Opx_par{index_Na2O(2)}*(1-DeltaNa2O);
X_Gt = X_Gt_par{index_Na2O(1)}*DeltaNa2O + X_Gt_par{index_Na2O(2)}*(1-DeltaNa2O);
X_Sp = X_Sp_par{index_Na2O(1)}*DeltaNa2O + X_Sp_par{index_Na2O(2)}*(1-DeltaNa2O);
X_Pl = X_Pl_par{index_Na2O(1)}*DeltaNa2O + X_Pl_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Melt = Vol_Melt_par{index_Na2O(1)}*DeltaNa2O + Vol_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Ol = Vol_Ol_par{index_Na2O(1)}*DeltaNa2O + Vol_Ol_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Cpx = Vol_Cpx_par{index_Na2O(1)}*DeltaNa2O + Vol_Cpx_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Opx = Vol_Opx_par{index_Na2O(1)}*DeltaNa2O + Vol_Opx_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Gt = Vol_Gt_par{index_Na2O(1)}*DeltaNa2O + Vol_Gt_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Sp = Vol_Sp_par{index_Na2O(1)}*DeltaNa2O + Vol_Sp_par{index_Na2O(2)}*(1-DeltaNa2O);
Vol_Pl = Vol_Pl_par{index_Na2O(1)}*DeltaNa2O + Vol_Pl_par{index_Na2O(2)}*(1-DeltaNa2O);

Rho_Melt = Rho_Melt_par{index_Na2O(1)}*DeltaNa2O + Rho_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
Rho_Ol = Rho_Ol_par{index_Na2O(1)}*DeltaNa2O + Rho_Ol_par{index_Na2O(2)}*(1-DeltaNa2O);
Rho_Cpx = Rho_Cpx_par{index_Na2O(1)}*DeltaNa2O + Rho_Cpx_par{index_Na2O(2)}*(1-DeltaNa2O);
Rho_Opx = Rho_Opx_par{index_Na2O(1)}*DeltaNa2O + Rho_Opx_par{index_Na2O(2)}*(1-DeltaNa2O);
Rho_Gt = Rho_Gt_par{index_Na2O(1)}*DeltaNa2O + Rho_Gt_par{index_Na2O(2)}*(1-DeltaNa2O);
Rho_Sp = Rho_Sp_par{index_Na2O(1)}*DeltaNa2O + Rho_Sp_par{index_Na2O(2)}*(1-DeltaNa2O);
Rho_Pl = Rho_Pl_par{index_Na2O(1)}*DeltaNa2O + Rho_Pl_par{index_Na2O(2)}*(1-DeltaNa2O);

H_Solid = H_Solid_par{index_Na2O(1)}*DeltaNa2O + H_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);
H_Melt = H_Melt_par{index_Na2O(1)}*DeltaNa2O + H_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
H_Gt = H_Gt_par{index_Na2O(1)}*DeltaNa2O + H_Gt_par{index_Na2O(2)}*(1-DeltaNa2O);
H_Sp = H_Sp_par{index_Na2O(1)}*DeltaNa2O + H_Sp_par{index_Na2O(2)}*(1-DeltaNa2O);

C_Solid = C_Solid_par{index_Na2O(1)}*DeltaNa2O + C_Solid_par{index_Na2O(2)}*(1-DeltaNa2O);
C_Melt = C_Melt_par{index_Na2O(1)}*DeltaNa2O + C_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
C_Cpx = C_Cpx_par{index_Na2O(1)}*DeltaNa2O + C_Cpx_par{index_Na2O(2)}*(1-DeltaNa2O);
C_Gt = C_Gt_par{index_Na2O(1)}*DeltaNa2O + C_Gt_par{index_Na2O(2)}*(1-DeltaNa2O);


if isnan(Cp_Solid)
    if isnan(Cp_Solid_par{index_Na2O(1)}) && isnan(Cp_Solid_par{index_Na2O(2)})
        Cp_Solid_vec = cell2mat(Cp_Solid_par);
        if any(~isnan(Cp_Solid_vec))
            index_no_nan = ~isnan(Cp_Solid_vec);
            Cp_Solid = sum(Cp_Solid_vec(index_no_nan))/sum(index_no_nan);
        end
    else
        if isnan(Cp_Solid_par{index_Na2O(1)})
            Cp_Solid = Cp_Solid_par{index_Na2O(2)};
        else
            Cp_Solid = Cp_Solid_par{index_Na2O(1)};
        end
    end
end
    
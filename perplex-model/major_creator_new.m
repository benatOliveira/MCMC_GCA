function [SiO2,Al2O3,FeOt,MnO,MgO,CaO,Na2O,Cr2O3,TiO2] = major_creator_new(Al2O3_ox,MnO_ox,Na2O_ox,Cr2O3_ox,TiO2_ox)
% This function uses the correlations explored and deduced in
% Afonso et al. (2013a); these are the same as those used in the LitMod
% suite of inversion codes.
% FeO is the most contrained a priori (i.e. most variable in reality)

% Al_ox needs to be inputed in wt% [valid range 0.1 < Al_ox < 4.5 ]
%--------------------------------------------------------------------------
Al2O3 = Al2O3_ox;
MnO = MnO_ox;
Na2O = Na2O_ox;
Cr2O3 = Cr2O3_ox;
TiO2 = TiO2_ox;

MgO = 49.369 + (-4.106 * Al2O3) + (0.343 * Al2O3.^2); %
CaO = -0.164 + (0.906 * Al2O3 );

FeOt(CaO > 1) = 7.4527 + (0.5689 * CaO(CaO > 1) ) + (-0.0863 * CaO(CaO > 1).^2);

FeOt(CaO <= 1) = 7.94 - (7.94-6.5) * (1-CaO(CaO <= 1))  ;

FeOt = FeOt';
SiO2 = 100 - Al2O3 - FeOt - MnO - MgO - CaO - Na2O - Cr2O3 - TiO2;

end
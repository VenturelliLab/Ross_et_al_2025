function dxdt = model_biomolecular(t, x, params, chemostat_params)
%% model_biomolecular
% 
% This function models the dynamics of two microbial populations growing 
% in a chemostat, where growth is coupled to glucose and cross-feeding of 
% two amino acids (e.g., tyrosine and phenylalanine).
%
% Inputs:
%   t - Current time (not used explicitly but required for ODE solvers)
%   x - State vector containing:
%           x(1) = N1, biomass concentration of strain 1
%           x(2) = N2, biomass concentration of strain 2
%           x(3) = R1, concentration of amino acid 1 (produced by N1, consumed by N2)
%           x(4) = R2, concentration of amino acid 2 (produced by N2, consumed by N1)
%           x(5) = R3, concentration of glucose
%   params - Vector of kinetic parameters:
%           [mu1_max, k13, k15, q13, q12, q11, k4, n4, 
%            mu2_max, k23, k24, q23, q21, q22, k5, n5]
%   chemostat_params - Vector of chemostat parameters:
%           [D, R1_in, R2_in, R3_in]
%           where D is the dilution rate, and R*_in are the influent concentrations
%
% Outputs:
%   dxdt - Time derivatives of the state variables, returned as a column vector

N1 = x(1); N2 = x(2); R1 = x(3); R2 = x(4); R3 = x(5);

mu1_max = params(01); k13 = params(02); k15 = params(03); q13 = params(04); q12 = params(05); q11 = params(06); k4 = params(07); n4 = params(08);
mu2_max = params(09); k23 = params(10); k24 = params(11); q23 = params(12); q21 = params(13); q22 = params(14); k5 = params(15); n5 = params(16);

D = chemostat_params(1);
R1_in = chemostat_params(2);
R2_in = chemostat_params(3);
R3_in = chemostat_params(4);

mu1 = mu1_max * R3/(k13 + R3) * R2/(k15 + R2);
mu2 = mu2_max * R3/(k23 + R3) * R1/(k24 + R1);


%% Equations
dN1dt = mu1*N1 - D*N1;
dN2dt = mu2*N2 - D*N2;
dR1dt = D*(R1_in - R1) + q11*(1/(1 + (R2/k4)^n4))*N1 - mu2*q21*N2;
dR2dt = D*(R2_in - R2) + q22*(1/(1 + (R1/k5)^n5))*N2 - mu1*q12*N1;
dR3dt = D*(R3_in - R3) - mu1*N1*q13 - mu2*N2*q23;

dxdt = [dN1dt;
        dN2dt;
        dR1dt;
        dR2dt;
        dR3dt];

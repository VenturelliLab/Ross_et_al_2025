function dxdt = model_reduced(t, x, p, c)
%% model_reduced
%
% This function models a simplified system describing the competition 
% between two microbial strains in a chemostat environment, with growth 
% coupled to glucose and two amino acids.
%
% Inputs:
%   t - Current time (required input for ODE solvers, though not explicitly used)
%   x - State vector containing:
%           x(1) = R1, concentration of amino acid 1
%           x(2) = R2, concentration of amino acid 2
%           x(3) = f, fraction of total biomass composed of strain 1
%   p - Parameter vector containing:
%           a12  - growth affinity coefficient for strain 1 with R2
%           a13  - growth affinity coefficient for strain 1 with R3
%           q12  - consumption coefficient of R2 by strain 1
%           q11  - production coefficient of R1 by strain 1
%           a21  - growth affinity coefficient for strain 2 with R1
%           a23  - growth affinity coefficient for strain 2 with R3
%           q21  - consumption coefficient of R1 by strain 2
%           q22  - production coefficient of R2 by strain 2
%           N_tot - total biomass concentration (constant)
%           R3    - concentration of glucose (assumed constant)
%   c - Chemostat parameter vector:
%           c(1) = D, dilution rate
%
% Outputs:
%   dxdt - Time derivatives of the state variables, returned as a column vector

R1 = x(1); R2 = x(2); f = x(3);


a12 = p(1); a13 = p(2); q12 = p(3); q11 = p(4);
a21 = p(5); a23 = p(6); q21 = p(7); q22 = p(8);

N_tot = p(9); R3 = p(10);

D = c(1);

mu12 = a12*R2; mu13 = a13*R3;
mu21 = a21*R1; mu23 = a23*R3;

mu1 = min([mu12, mu13]);
mu2 = min([mu21, mu23]);


%% Equations
dR1dt = -D*R1 + q11*(mu13 - mu1)*f*N_tot - q21*mu2*(1-f)*N_tot;
dR2dt = -D*R2 + q22*(mu23 - mu2)*(1-f)*N_tot - q12*mu1*f*N_tot;
dfdt  = (mu1 - mu2)*f*(1-f);

dxdt = [dR1dt;
        dR2dt;
        dfdt];

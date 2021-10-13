function [AMatrix,BMatrix,DMatrix,SMatrix,Inertia] = ...
    reddyFgmMaterial(thickness,n)

%%% REDDY FGM MATERIAL

% plate thickness
h = thickness;

% elastic moduli
E1 = 14.4e9; E2 = 1.44e9;
rho1 = 12.2e3; rho2 = 1.22e3;
poisson = 0.38;
kapa = 5*(1+poisson)/(6+5*poisson);

c1 = 1/(1-poisson^2);
c2 = poisson/(1-poisson^2);
c3 = 1/2/(1+poisson);

% stiffness calculation
A11 = c1*((E1-E2)*h/(n+1) + h*E2);
A12 = c2*((E1-E2)*h/(n+1) + h*E2);
A66 = c3*((E1-E2)*h/(n+1) + h*E2);

A44 = kapa*c3*((E1-E2)*h/(n+1) + h*E2);
A55 = A44;

B11 = c1*(E1-E2)*n*h^2/(2*(n+1)*(n+2));
B12 = c2*(E1-E2)*n*h^2/(2*(n+1)*(n+2));
B66 = c3*(E1-E2)*n*h^2/(2*(n+1)*(n+2));

D11 = c1*((E1-E2)*h^3*(2+n+n^2)/(4*(n+1)*(n+2)*(n+3)) + E2*h^3/12);
D12 = c2*((E1-E2)*h^3*(2+n+n^2)/(4*(n+1)*(n+2)*(n+3)) + E2*h^3/12);
D66 = c3*((E1-E2)*h^3*(2+n+n^2)/(4*(n+1)*(n+2)*(n+3)) + E2*h^3/12);

% inertia calculation
I0 = (rho1-rho2)*h/(n+1) + rho2*h;
I1 = (rho1-rho2)*n*h^2/(2*(n+1)*(n+2));
I2 = (rho1-rho2)*h^3*(2+n+n^2)/(4*(n+1)*(n+2)*(n+3)) + rho2*h^3/12;


AMatrix = [A11,A12,0;A12,A11,0;0,0,A66];
BMatrix = [B11,B12,0;B12,B11,0;0,0,B66];
DMatrix = [D11,D12,0;D12,D11,0;0,0,D66];
SMatrix = [A44,0;0,A55];

Inertia = [I0 0 0 0 0;
           0 I2 0 I1 0;
           0 0 I2 0 I1;
           0 I1 0 I0 0;
           0 0 I1 0 I0];
end
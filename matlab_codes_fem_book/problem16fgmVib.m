% ................................................................
% MATLAB codes for Finite Element Analysis
% problem16fgmVib.m
% Functionally graded Timoshenko beam in free vibrations 
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% E1: modulus of elasticity of material 1
% E2: modulus of elasticity of material 2
% rho1: density of material 1
% rho2: density of material 2
% L: length of beam
% thickness: height of cross-section
% width: width of the cross-section
E1 = 14.4e9; E2 = 1.44e9;
rho1 = 12.2e3; rho2 = 1.22e3;
poisson = 0.38;
thickness = 17.6e-6;
width = 2*thickness;

L = 20*thickness;

n = 10; % FGM power-law index

M = E1/E2;
A0 = width*thickness;
B0 = width*thickness^2;
I0 = width*thickness^3/12;

Axx = E2*A0*(M+n)/(1+n);
Bxx = E2*B0*n*(M-1)/2/(1+n)/(2+n);
Dxx = E2*I0*((6+3*n+3*n^2)*M + 8*n+3*n^2+n^3)/(6+11*n+6*n^2+n^3);

kapa = 5*(1+poisson)/(6+5*poisson);
Sxz = kapa*E2*A0/2/(1+poisson)*(M+n)/(n+1);

% constitutive matrix
C = [Axx Bxx 0; Bxx Dxx 0; 0 0 Sxz];

m0 = A0*(rho1 + n*rho2)/(n+1);
m1 = B0*n*(rho1-rho2)/2/(n+1)/(n+2);
m2 = I0*((6+3*n+3*n^2)*rho1 + (8*n+3*n^2+n^3)*rho2)/ ...
    (6+11*n+6*n^2+n^3);
% inertia matrix
I = [m0 0 m1; 0 m0 0; m1 0 m2];

% mesh
numberElements = 40;  
nodeCoordinates = linspace(0,L,numberElements+1);
xx = nodeCoordinates';
elementNodes = zeros(size(nodeCoordinates,2)-1,2);
for i = 1:size(nodeCoordinates,2)-1
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1;
end
% generation of coordinates and connectivities
numberNodes = size(xx,1);

% GDof: global number of degrees of freedom
GDof = 3*numberNodes; 

% computation of the system stiffness matrix 
[stiffness,~,mass] = ...
    formStiffnessMassTimoshenkoFgmBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,C,1,I,thickness);

% boundary conditions (simply-supported at both bords)
fixedNodeU = [1];
fixedNodeW = [1 ; numberNodes];
fixedNodeTX = [];
prescribedDof = [fixedNodeU; fixedNodeW+numberNodes; ...
    fixedNodeTX+2*numberNodes];

% free vibrations
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);

omega = sqrt(eigenvalues)*L*L*sqrt(rho2*A0/E2/I0);
% display first 2 dimensionless frequencies
omega(1:3)

% drawing mesh and deformed shape
modeNumber = 4;
V1 = modes(:,1:modeNumber);

% drawing eigenmodes
figure
drawEigenmodes1D(modeNumber,numberNodes,...
    V1(1+numberNodes:2*numberNodes,:),xx)
% ................................................................
% MATLAB codes for Finite Element Analysis
% problem16vibrations.m
% Timoshenko beam in free vibrations
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E: modulus of elasticity
% G: shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
E = 1e8; poisson = 0.30; L = 1; thickness = 0.001; rho = 1;
I = thickness^3/12; EI = E*I; A = 1*thickness; kapa = 5/6;

% constitutive matrix
G = E/2/(1+poisson);
C = [EI 0; 0  kapa*thickness*G];

% mesh
numberElements = 50;
nodeCoordinates = linspace(0,L,numberElements+1);
xx = nodeCoordinates'; x = xx';
elementNodes = zeros(size(nodeCoordinates,2)-1,2);
for i = 1:size(nodeCoordinates,2)-1
    elementNodes(i,1) = i; 
    elementNodes(i,2) = i+1;
end
% generation of coordinates and connectivities
numberNodes = size(xx,1);

% GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% computation of the system stiffness, force, mass
[stiffness,force,mass] = ...
    formStiffnessMassTimoshenkoBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,C,0,rho,I,thickness);

% boundary conditions (simply-supported at both ends)
% fixedNodeW  = [1 ; numberNodes];
% fixedNodeTX = []; 
% boundary conditions (clamped at both ends)
% fixedNodeW  = [1 ; numberNodes];
% fixedNodeTX = fixedNodeW; 
% boundary conditions (cantilever)
fixedNodeW  = [1];
fixedNodeTX = [1];
prescribedDof = [fixedNodeW; fixedNodeTX+numberNodes];

% free vibration problem
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);

omega = sqrt(eigenvalues)*L*L*sqrt(rho*A/E/I);
% display first 2 dimensionless frequencies
omega(1:2)

% drawing mesh and deformed shape
modeNumber = 4;
V1 = modes(:,1:modeNumber);

% drawing eigenmodes
figure
drawEigenmodes1D(modeNumber,numberNodes,V1,x)
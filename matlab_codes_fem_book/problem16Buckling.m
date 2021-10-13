% ................................................................
% MATLAB codes for Finite Element Analysis
% problem16Buckling.m
% Timoshenko beam under buckling loads (P)
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E: modulus of elasticity
% G: shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
E = 10e6; poisson = 0.333; L = 1; thickness = 0.1;
I = thickness^3/12; EI = E*I; A = 1*thickness; kapa = 5/6;

% constitutive matrix
G = E/2/(1+poisson);
C = [EI 0; 0  kapa*thickness*G];

% mesh
numberElements = 40;
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

% computation of the system stiffness, Kg
[stiffness,Kg] = ...
    formStiffnessBucklingTimoshenkoBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,C,I,thickness);

% boundary conditions (CC)
fixedNodeW = find(xx==min(nodeCoordinates(:)) ...
    | xx==max(nodeCoordinates(:)));
fixedNodeTX = fixedNodeW;
prescribedDof = [fixedNodeW; fixedNodeTX+numberNodes];

% boundary conditions (SS)
% fixedNodeW = find(xx==min(nodeCoordinates(:))...
%     | xx==max(nodeCoordinates(:)));
% prescribedDof = [fixedNodeW];

% buckling problem
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
    stiffness,Kg,0);

% reordering eigenvalues
[eigenvalues,ii] = sort(eigenvalues);  modes = modes(:,ii);

% Bazant & Cedolin solution for SS and CC
PcrSS = pi*pi*E*I/L^2*(1/(1+pi*pi*E*I/(L*L*kapa*G*A)))
PcrCC = pi*pi*E*I/(L/2)^2*(1/(1+pi*pi*E*I/(L*L/4*kapa*G*A)))

% drawing mesh and deformed shape
modeNumber = 4;
V1 = modes(:,1:modeNumber);

% drawing eigenmodes
figure
drawEigenmodes1D(modeNumber,numberNodes,V1,x)
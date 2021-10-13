% ................................................................
% MATLAB codes for Finite Element Analysis
% problem16.m
% Timoshenko beam in bending
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E: modulus of elasticity
% G: shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
E = 1e8; poisson = 0.30; L = 1; thickness = 0.001;
I = thickness^3/12;
EI = E*I;
kapa = 5/6;

P = -1; % uniform pressure
% constitutive matrix
G = E/2/(1+poisson);
C = [EI 0; 0 kapa*thickness*G];

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
GDof = 2*numberNodes; 

% stiffness matrix and force vector
[stiffness,force] = ...
    formStiffnessMassTimoshenkoBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,C,P,1,I,thickness);

% boundary conditions (simply-supported at both ends)
% fixedNodeW = [1 ; numberNodes];
% fixedNodeTX = []; 
% boundary conditions (clamped at both ends)
% fixedNodeW = [1 ; numberNodes];
% fixedNodeTX = fixedNodeW; 
% boundary conditions (cantilever)
fixedNodeW = [1];
fixedNodeTX = [1];
prescribedDof = [fixedNodeW; fixedNodeTX+numberNodes];

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    GDof,prescribedDof)

U = displacements;
ws = 1:numberNodes;

% max displacement
disp('max displacement')
min(U(ws))


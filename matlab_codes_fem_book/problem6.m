%................................................................
% MATLAB codes for Finite Element Analysis
% problem6.m
% ref: D. Logan, A first course in the finite element method,
% third Edition, mixing trusses with springs
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% A: area of cross section
E = 210000; A = 500; EA = E*A;

% generation of coordinates and connectivities
nodeCoordinates = [0 0;-5000*cos(pi/4) 5000*sin(pi/4); -10000 0];
elementNodes = [1 2;1 3;1 4];
numberElements = size(elementNodes,1);
numberNodes = size(nodeCoordinates,1)+1; % spring added
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

% for structure:
%   displacements: displacement vector
%   force : force vector
%   stiffness: stiffness matrix
GDof = 2*numberNodes;
U = zeros(GDof,1);
force = zeros(GDof,1);
stiffness = zeros(GDof);

% applied load at node 2
force(2) = -25000;

% computation of the system stiffness matrix
[stiffness] = ...
    formStiffness2Dtruss(GDof,numberElements-1, ...
    elementNodes,numberNodes,nodeCoordinates,xx,yy,EA);

% spring stiffness in global Dof
stiffness([2 7],[2 7]) = stiffness([2 7],[2 7]) + 2000*[1 -1;-1 1];

% boundary conditions and solution
prescribedDof = (3:8)';

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    GDof,prescribedDof)

% stresses at elements
stresses2Dtruss(numberElements-1,elementNodes, ...
    xx,yy,displacements,E)


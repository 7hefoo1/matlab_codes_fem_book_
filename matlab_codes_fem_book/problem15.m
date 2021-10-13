% ................................................................
% MATLAB codes for Finite Element Analysis
% problem15.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E: modulus of elasticity
% I: second moments of area
% J: polar moment of inertia
% G: shear modulus
E = 210e9; G = 84e9; I = 20e-5; J = 5e-5;

% generation of coordinates and connectivities
nodeCoordinates = [4 4; 0 4; 0 0 ; 4 0];
xx = nodeCoordinates(:,1);   
yy = nodeCoordinates(:,2); 
elementNodes = [1 2; 3 1; 4 1]; 
numberNodes = size(nodeCoordinates,1);
numberElements = size(elementNodes,1);

% GDof: global number of degrees of freedom
GDof = 3*numberNodes; 

%force vector
force = zeros(GDof,1);
force(1) = -20e3;

% stiffness matrix
stiffness = formStiffnessGrid(GDof,numberElements, ...
    elementNodes,xx,yy,E,I,G,J);

% boundary conditions
prescribedDof = (4:12)';

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    GDof,prescribedDof)

% forces in elements
disp('forces in elements')
EF = forcesInElementGrid(numberElements,elementNodes, ...
    xx,yy,E,I,G,J,displacements)
% ................................................................
% MATLAB codes for Finite Element Analysis
% problem13vib.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% I: second moments of area
% J: polar moment of inertia
% G: shear modulus
E = 210e9; A = 0.02; rho = 7850;
Iy = 10e-5; Iz = 20e-5; J = 5e-5; G = 84e9;

% generation of coordinates and connectivities
nodeCoordinates = [0 0 0; 
    0 0 4; 
    4 0 4; 
    4 0 0;
    0 5 0; 
    0 5 4; 
    4 5 4; 
    4 5 0];
xx = nodeCoordinates(:,1);   
yy = nodeCoordinates(:,2); 
zz = nodeCoordinates(:,3);
elementNodes = [1 5;2 6;3 7; 4 8; 5 6; 6 7; 7 8; 8 5]; 
numberNodes = size(nodeCoordinates,1);
numberElements = size(elementNodes,1);

% for structure:
%   displacements: displacement vector
%   stiffness: stiffness matrix
%   mass: mass matrix
%   GDof: global number of degrees of freedom
GDof = 6*numberNodes; 
U = zeros(GDof,1);

% stiffness matrix
[stiffness] = ...
    formStiffness3Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);

% mass matrix 
[mass] = ...
    formMass3Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,rho,A,Iz,Iy);

% boundary conditions
prescribedDof = 1:24;

% solution
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);

omega = sqrt(eigenvalues);

% drawing mesh and deformed shape
modeNumber = 1;
U = modes(:,modeNumber);

figure
XX = U(1:6:6*numberNodes);
YY = U(2:6:6*numberNodes);
ZZ = U(3:6:6*numberNodes);
scaleFact = 20;
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY ZZ],elementNodes, ...
    'L2','k.');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');
axis equal
set(gca,'fontsize',18)
view(170,-45)

% plot interpolated deformed shape acoording 
% to Lagrange and Hermite shape functions
drawInterpolatedFrame3D
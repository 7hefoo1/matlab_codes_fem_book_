%% ................................................................
% MATLAB codes for Finite Element Analysis
% problem11bvib.m
% 2D frame
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear
close all

% E; modulus of elasticity
% I: second moment of area
E = 210000; A = 200; I = 2e8; rho = 8.05e-9; 
EA = E*A; EI = E*I; rhoA = rho*A;

% generation of coordinates and connectivities
numberElements = 12;
nodeCoordinates = [0 0;0 1500;0 3000;0 4500 ; 
                   0 6000;1500 6000;3000 6000;
                   4500 6000;6000 6000;6000 4500;
                   6000 3000;6000 1500;6000 0];
elementNodes = zeros(numberElements,2);
for i = 1:numberElements
    elementNodes(i,1) = i;
    elementNodes(i,2) = i+1;
end
numberNodes = size(nodeCoordinates,1);
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

% for structure:
%   displacements: displacement vector
%   stiffness: stiffness matrix
%   GDof: global number of degrees of freedom
GDof = 3*numberNodes; 
U = zeros(GDof,1);

% stiffness matrix
[stiffness] = ...
    formStiffness2Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,xx,yy,EI,EA);

% mass matrix
[mass] = ...
    formMass2Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,xx,yy,rhoA);

% boundary conditions and solution
prescribedDof = [1 13 14 26 27 39]';

% solution
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);

omega = sqrt(eigenvalues);

% drawing mesh and deformed shape
modeNumber = 3;
U = modes(:,modeNumber);

figure
XX = U(1:numberNodes); YY = U(numberNodes+1:2*numberNodes);
dispNorm = max(sqrt(XX.^2+YY.^2));
scaleFact = 20*dispNorm;
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY],elementNodes, ...
    'L2','k.');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');
axis equal
set(gca,'fontsize',18)

% plot interpolated deformed shape acoording 
% to Lagrange and Hermite shape functions
displacements = U;
drawInterpolatedFrame2D
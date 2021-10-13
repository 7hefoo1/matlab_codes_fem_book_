%................................................................
% MATLAB codes for Finite Element Analysis
% problem5vib.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% A: area of cross section
E = 70000; A = 300; EA = E*A; rho = 1000; rhoA = rho*A;

% generation of coordinates and connectivities
elementNodes = [1 2;1 3;2 3;2 4;1 4;3 4;3 6;4 5;4 6;3 5;5 6];
nodeCoordinates = [0 0;0 3000;3000 0;3000 3000;6000 0;6000 3000];
numberElements = size(elementNodes,1);
numberNodes = size(nodeCoordinates,1);
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

% for structure:
%   displacements: displacement vector
%   force : force vector
%   stiffness: stiffness matrix
GDof = 2*numberNodes;
U = zeros(GDof,1);

% computation of the system stiffness matrix
[stiffness] = ...
    formStiffness2Dtruss(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,xx,yy,EA);

% computation of the system stiffness matrix
[mass] = ...
    formMass2Dtruss(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,xx,yy,rhoA);

% boundary conditions and solution
prescribedDof = [1 2 10]';

% free vibration problem
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);
us = 1:2:2*numberNodes-1;
vs = 2:2:2*numberNodes;

modeNumber = 1;

% drawing displacements
figure
L = xx(2)-xx(1);
XX = modes(us,modeNumber); YY = modes(vs,modeNumber);
dispNorm = max(sqrt(XX.^2+YY.^2));
scaleFact = 1e12*dispNorm;
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY], ...
    elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');
axis equal
set(gca,'fontsize',18)

omega = sqrt(eigenvalues)
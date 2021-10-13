% ................................................................
% MATLAB codes for Finite Element Analysis
% problem7vib.m
% A 3D truss example in free vibrations
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% A: area of cross section
E = 1.2e6; 
A = [0.302;0.729;0.187]; % area for various sections
rho = [1;1;1]; % density for various sections

% generation of coordinates and connectivities
nodeCoordinates = [72 0 0; 0 36 0;  0 36 72; 0 0 -48];
elementNodes = [1 2;1 3;1 4]; 
numberElements = size(elementNodes,1);
numberNodes = size(nodeCoordinates,1); 
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

% for structure:
%   displacements: displacement vector
%   force : force vector
%   stiffness: stiffness matrix
%   GDof: global number of degrees of freedom
GDof = 3*numberNodes; 
U = zeros(GDof,1);
force = zeros(GDof,1);

% stiffness matrix
[stiffness] = ...
    formStiffness3Dtruss(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,E,A);

% mass matrix
[mass] = ...
    formMass3Dtruss(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,rho,A);

% boundary conditions and solution
prescribedDof = [2 4:12]';

% free vibration problem
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);
us = 1:3:3*numberNodes-2;
vs = 2:3:3*numberNodes-1;
ws = 3:3:3*numberNodes;

modeNumber = 1;

% drawing displacements
figure
L = xx(2)-xx(1);
XX = modes(us,modeNumber);
YY = modes(vs,modeNumber);
ZZ = modes(ws,modeNumber);
dispNorm = max(sqrt(XX.^2+YY.^2+ZZ.^2));
scaleFact = 1e3*dispNorm;
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY ZZ], ...
    elementNodes,'L3','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L3','k.--');
axis equal
set(gca,'fontsize',18)
view(45,45)

omega = sqrt(eigenvalues)
% ................................................................
% MATLAB codes for Finite Element Analysis
% problem11b.m
% 2D frame
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% I: second moment of area
E = 210000; A = 200; I = 2e8; EA = E*A; EI = E*I;

% generation of coordinates and connectivities
numberElements = 12;
nodeCoordinates = [0 0;0 1500;0 3000;0 4500;
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
%   force: force vector
%   stiffness: stiffness matrix
%   GDof: global number of degrees of freedom
GDof = 3*numberNodes; 
force = zeros(GDof,1);
stiffness = zeros(GDof); 

%force vector
force(5) = 15000;
force(31) = 10e6;

% stiffness matrix
[stiffness] = ...
    formStiffness2Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,xx,yy,EI,EA);

% boundary conditions and solution
prescribedDof = [1 13 14 26 27 39]';

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    GDof,prescribedDof)

%drawing mesh and deformed shape
figure
XX = displacements(1:numberNodes);
YY = displacements(numberNodes+1:2*numberNodes);
dispNorm = max(sqrt(XX.^2+YY.^2));
scaleFact = 50*dispNorm;
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY],elementNodes, ...
    'L2','k.');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');
axis equal
set(gca,'fontsize',18)

% plot interpolated deformed shape acoording 
% to Lagrange and Hermite shape functions
drawInterpolatedFrame2D
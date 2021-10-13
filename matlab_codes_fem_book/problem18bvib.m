% ................................................................
% MATLAB codes for Finite Element Analysis
% problem18vib.m
% 2D problem: beam in free vibrations using Q9 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
E = 10e7; poisson = 0.30; rho = 1000;

% matrix C
C = E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

% mesh generation
Lx = 5;
Ly = 1;
numberElementsX = 20;
numberElementsY = 10;
numberElements = numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(Lx,Ly,numberElementsX,numberElementsY,'Q9');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q9','-');
axis equal

numberNodes = size(xx,1);
% GDof: global number of degrees of freedom
GDof = 2*numberNodes;

% stiffness and mass matrices
[stiffness,mass] = formStiffnessMass2D(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,C,rho,1,...
    'Q9','complete');

% boundary conditions
fixedNodeX = find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY = find(nodeCoordinates(:,1)==0);  % fixed in YY
prescribedDof = [fixedNodeX; fixedNodeY+numberNodes];


% solution
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
    stiffness,mass,15);

omega = sqrt(eigenvalues);
% sort out eigenvalues
[omega,ii] = sort(omega);
modes = modes(:,ii);

% drawing mesh and deformed shape
modeNumber = 3;
displacements = modes(:,modeNumber);

% displacements and deformed shape
UX = displacements(1:numberNodes);
UY = displacements(numberNodes+1:GDof);
scaleFactor = 0.5;

% deformed shape
figure
drawingField(nodeCoordinates+scaleFactor*[UX UY], ...
    elementNodes,'Q9',UX);%U XX
hold on
drawingMesh(nodeCoordinates+scaleFactor*[UX UY], ...
    elementNodes,'Q9','-');
drawingMesh(nodeCoordinates,elementNodes,'Q9','--');
colorbar
title('Displacement field u_x (on deformed shape)')
axis off
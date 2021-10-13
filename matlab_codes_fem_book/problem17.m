% ................................................................
% MATLAB codes for Finite Element Analysis
% problem17.m
% 2D problem: thin plate in tension using Q4 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% material properties
E = 10e7; poisson = 0.30;  

% matrix C
C = E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

% load
P = 1e6;

% mesh generation
Lx = 5;
Ly = 1;
numberElementsX = 10;
numberElementsY = 5;
numberElements = numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(Lx,Ly,numberElementsX,numberElementsY,'Q4');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q4','-');
axis equal

numberNodes = size(xx,1);
% GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% calculation of the system stiffness matrix
stiffness = formStiffnessMass2D(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,C,1,1,'Q4','complete');

% boundary conditions 
fixedNodeX = find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY = find(nodeCoordinates(:,2)==0);  % fixed in YY
prescribedDof = [fixedNodeX; fixedNodeY+numberNodes];

% force vector (distributed load applied at xx=Lx)
force = zeros(GDof,1);
rightBord = find(nodeCoordinates(:,1)==Lx);
force(rightBord) = P*Ly/numberElementsY;
force(rightBord(1)) = P*Ly/numberElementsY/2;
force(rightBord(end)) = P*Ly/numberElementsY/2;

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% displacements
disp('Displacements')
jj = 1:GDof; format
f = [jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)
UX = displacements(1:numberNodes);
UY = displacements(numberNodes+1:GDof);
scaleFactor = 10;

% deformed shape
figure
drawingField(nodeCoordinates+scaleFactor*[UX UY], ...
    elementNodes,'Q4',UX);%U XX
hold on
drawingMesh(nodeCoordinates+scaleFactor*[UX UY], ...
    elementNodes,'Q4','-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','--');
colorbar
title('Displacement field u_x (on deformed shape)')
axis off

% stresses at nodes
[stress,strain] = stresses2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,displacements,...
    C,'Q4','complete');

% drawing stress fields on deformed shape
figure; hold on;
drawingField2(nodeCoordinates,elementNodes,...
    scaleFactor,UX,UY,'Q4',stress(:,:,1))
axis equal
drawingMesh(nodeCoordinates,elementNodes,'Q4','--');
colorbar
title('Stress field \sigma_{xx} (on deformed shape)')
axis off
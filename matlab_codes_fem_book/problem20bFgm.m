% ................................................................
% MATLAB codes for Finite Element Analysis
% problem20Fgm.m
% functionally graded plate using Q9 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
thickness = 17.6e-6;
n = 10; % power-law index

% load
P = -1;

% mesh generation
L = 20*thickness;
numberElementsX = 10;
numberElementsY = 10;
numberElements = numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY,'Q9');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q9','-');
axis equal

numberNodes = size(xx,1);

% GDof: global number of degrees of freedom
GDof = 5*numberNodes;

% computation of the system stiffness matrix
% the shear correction factors are automatically
% calculted for any laminate
[AMatrix,BMatrix,DMatrix,SMatrix] = reddyFgmMaterial(thickness,n);

stiffness = formStiffnessMatrixMindlinlaminated5dof ...
    (GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,AMatrix, ...
    BMatrix,DMatrix,SMatrix,'Q9','third','complete');

% computation of the system force vector
[force] = ...
    formForceVectorMindlin5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P,'Q9','complete');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
U = solution(GDof,prescribedDof,stiffness,force);

% drawing deformed shape and normalize results
% to compare with Srinivas
ws = 1:numberNodes;
disp('maximum displacement')
abs(min(U(ws))*1.44e9*thickness^3/P/L^4)

% surface representation
figure; hold on
for k = 1:size(elementNodes,1)
    patch(nodeCoordinates(elementNodes(k,1:4),1),...
        nodeCoordinates(elementNodes(k,1:4),2),...
        U(elementNodes(k,1:4)),...
        U(elementNodes(k,1:4)))
end
set(gca,'fontsize',18)
view(45,45)
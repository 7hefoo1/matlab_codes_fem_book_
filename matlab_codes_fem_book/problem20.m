% ................................................................
% MATLAB codes for Finite Element Analysis
% problem20.m
% laminated plate: Srinivas problem using Q4 elements
% S. Srinivas, A refined analysis of composite laminates,
% J. Sound and Vibration, 30 (1973),495--507.
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
thickness = 0.1;

% load
P = -1;

% mesh generation
L = 1;
numberElementsX = 10;
numberElementsY = 10;
numberElements = numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY,'Q4');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q4','-');
axis equal

numberNodes = size(xx,1);

% GDof: global number of degrees of freedom
GDof = 5*numberNodes;

% computation of the system stiffness matrix
% the shear correction factors are automatically
% calculted for any laminate
[AMatrix,BMatrix,DMatrix,SMatrix,qbarra] = ...
    srinivasMaterial(thickness);

stiffness = formStiffnessMatrixMindlinlaminated5dof ...
    (GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,AMatrix, ...
    BMatrix,DMatrix,SMatrix,'Q4','complete','reduced');

% computation of the system force vector
[force] = ...
    formForceVectorMindlin5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P,'Q4','reduced');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
U = solution(GDof,prescribedDof,stiffness,force);

% drawing deformed shape and normalize results
% to compare with Srinivas
ws = 1:numberNodes;
disp('maximum displacement')
abs(min(U(ws))*0.999781/thickness)

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

% stress computation (Srinivas only)
disp('stress computation (Srinivas only)')
[stress_layer1,stress_layer2,...
    stress_layer3,shear_layer1,...
    shear_layer2] = SrinivasStress(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,...
    qbarra,U,thickness,'Q4','complete','reduced');


% normalized stresses, look for table in the book
format
[ abs(min(stress_layer3(:,3,1))),...
    abs(min(stress_layer2(:,3,1))), ...
    abs(min(stress_layer1(:,3,1))),...
    max(shear_layer2(:,1,1)),...
    max(shear_layer1(:,1,1))]
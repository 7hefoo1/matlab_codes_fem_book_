% ................................................................
% MATLAB codes for Finite Element Analysis
% problem20bBuckling.m
% buckling laminated plate using Q9 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
thickness = 0.001;

% initial stress matrix
sigmaX = 1/thickness;
sigmaXY = 0;
sigmaY = 0;
sigmaMatrix = [sigmaX sigmaXY; sigmaXY sigmaY];

%Mesh generation
Lx = 1;
Ly = 1;
numberElementsX = 10;
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
GDof = 5*numberNodes;

% computation of the system stiffness matrix
% the shear correction factors are automatically
% calculted for any laminate
[AMatrix,BMatrix,DMatrix,SMatrix] = reddyLaminateMaterialBuk(thickness);

stiffness = formStiffnessMatrixMindlinlaminated5dof ...
    (GDof,numberElements,elementNodes,numberNodes,nodeCoordinates, ...
    AMatrix, BMatrix,DMatrix,SMatrix,'Q9','third','complete');

% computation of the system force vector
geometric = formGeometricStiffnessMindlinlaminated5dof ...
    (GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,...
    sigmaMatrix,thickness,'Q9','complete','complete');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
% buckling analysis ...
[modes,lambda] = eigenvalue(GDof,prescribedDof,...
    stiffness,geometric,15);


% sort out eigenvalues
[lambda,ii] = sort(lambda);
modes = modes(:,ii);

% dimensionless omega
lambda_bar = lambda*Ly^2/pi/pi/DMatrix(2,2);
% lambda_bar = lambda*Ly^2/thickness^3
lambda_bar(1)

%drawing mesh and deformed shape
modeNumber = 1;
displacements = modes(:,modeNumber);

% surface representation
figure; hold on
for k = 1:size(elementNodes,1)
    patch(nodeCoordinates(elementNodes(k,1:4),1),...
        nodeCoordinates(elementNodes(k,1:4),2),...
        displacements(elementNodes(k,1:4)),...
        displacements(elementNodes(k,1:4)))
end
set(gca,'fontsize',18)
view(45,45)
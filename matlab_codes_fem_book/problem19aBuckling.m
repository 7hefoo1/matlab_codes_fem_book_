% ................................................................
% MATLAB codes for Finite Element Analysis
% problem19aBuckling.m
% Buckling analysis of Mindlin plates using Q8 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; colordef white

% material
E = 10920; poisson = 0.30; kapa = 5/6;
thickness = 0.001;
I = thickness^3/12;

% L: side lenght
L = 1;

% constitutive matrix
% bending part
C_bending = I*E/(1-poisson^2)*[1       poisson 0;
                               poisson 1       0;
                               0       0       (1-poisson)/2];
% shear part
C_shear = kapa*thickness*E/2/(1+poisson)*eye(2);

% initial stress matrix
sigmaX = 1/thickness;
sigmaXY = 0;
sigmaY = 0;
sigmaMatrix = [sigmaX sigmaXY; sigmaXY sigmaY];

% mesh generation ...
% numberElementsX: number of elements in x
% numberElementsY: number of elements in y
numberElementsX = 10;
numberElementsY = 10;
% number of elements
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY,'Q8');
xx = nodeCoordinates(:,1);   yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q8','-');
axis equal

numberNodes = size(xx,1);    % number of nodes
GDof = 3*numberNodes;        % total number of DOFs

% stiffness and geometric stiffness matrices
[stiffness] = ...
    formStiffnessMatrixMindlin(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,C_shear, ...
    C_bending,'Q8','third','complete');

[geometric] = ...
    formGeometricStiffnessMindlin(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,...
    thickness,'Q8','complete','complete');

% Essential boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% buckling analysis ...
[modes,lambda] = eigenvalue(GDof,prescribedDof,...
    stiffness,geometric,15);

% sort out eigenvalues
[lambda,ii] = sort(lambda);
modes = modes(:,ii);

% dimensionless omega
lambda_bar = lambda*L*L/pi/pi/C_bending(1,1)

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
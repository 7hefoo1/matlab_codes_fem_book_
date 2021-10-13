% ................................................................
% MATLAB codes for Finite Element Analysis
% problem20timeReddy.m
% laminated plate time transient using Q4 elements
% ref: J.N. Reddy, Mechanics of Laminated Composite Plates and 
% Shells 2nd Ed.
% Section 6.7.4 page 364
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
thickness = 2.5;

%Mesh generation
L = 25;
numberElementsX = 10; numberElementsY = 10;
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
[AMatrix,BMatrix,DMatrix,SMatrix,Inertia] = ...
    reddyLaminateMaterial(thickness);

stiffness = formStiffnessMatrixMindlinlaminated5dof ...
    (GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,AMatrix, ...
    BMatrix,DMatrix,SMatrix,'Q4','complete','reduced');

% computation of the system force vector
[force] = ...
    formForceVectorMindlin5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,0,'Q4','complete');

[mass] = ...
    formMassMatrixFgmPlate5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,Inertia,...
    'Q4','complete');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% Time transient simulation
timeStep = 5;
totalTime = 1000;
time = timeStep:timeStep:totalTime;

% Newton's parameters
alpha = 1/2; gamma = 1/2;

a1 = (1-alpha)*timeStep;
a2 = alpha*timeStep;
a3 = 2/(gamma*timeStep^2);
a4 = a3*timeStep;
a5 = (1-gamma)/gamma;

% initialization
displacementsTime = zeros(GDof,length(time));
velocitiesTime = zeros(GDof,length(time));
accelerationsTime = zeros(GDof,length(time));

% initial conditions
accelerationsTime(:,1) = mass\(force ...
    - stiffness*displacementsTime(:,1));

q0 = -1;
P = q0.*(1 - cos(0.0185*time));

for i = 2:length(time)
    [force] = ...
    formForceVectorMindlin5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P(i),'Q4','complete');
    
    forceHat = force + mass*(a3.*displacementsTime(:,i-1) ...
        + a4.*velocitiesTime(:,i-1) ...
        + a5.*accelerationsTime(:,i-1));
    stiffnessHat = stiffness + a3.*mass;
    
    displacementsTime(:,i) = solution(GDof,prescribedDof,...
        stiffnessHat,forceHat);
    
    accelerationsTime(:,i) = a3*(displacementsTime(:,i) ...
        - displacementsTime(:,i-1)) ...
        - a4.*velocitiesTime(:,i-1) ...
        - a5.*accelerationsTime(:,i-1);
    velocitiesTime(:,i) = velocitiesTime(:,i-1) ...
        + a1.*accelerationsTime(:,i-1) ...
        + a2.*accelerationsTime(:,i);
end

% dimensionless transverse displacement vs time
centralPt = find(xx==L/2 & yy==L/2);
w_bar = displacementsTime(centralPt,:)*2.1e6*thickness^3/q0/L^4*1e2;

figure;
hold on
plot(time,w_bar,'.-','linewidth',2,'markersize',16)
xlabel('time'); ylabel('central point motion')
ylim([-0.2 4.2])
set(gca,'linewidth',2,'fontsize',14)
grid on; box on;
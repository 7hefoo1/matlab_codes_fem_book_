% ................................................................
% MATLAB codes for Finite Element Analysis
% problem16timeReddy.m
% Timoshenko beam time transient analysis
% ref: J.N. Reddy, an introduction to Finite Element Method 3rd Ed.
% Example 6.2.2 page 332
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E: modulus of elasticity
% G: shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
poisson = 0.25; L = 1; thickness = 0.01; I = thickness^3/12;
E = 1/I; rho = 100; EI = E*I; kapa = 5/6;

% constitutive matrix
G = E/2/(1+poisson);
C = [EI 0; 0 kapa*thickness*G];

% mesh
numberElements = 40;
nodeCoordinates = linspace(0,L,numberElements+1);
xx = nodeCoordinates';
elementNodes = zeros(size(nodeCoordinates,2)-1,2);
for i = 1:size(nodeCoordinates,2)-1
    elementNodes(i,1)=i;
    elementNodes(i,2)=i+1;
end
% generation of coordinates and connectivities
numberNodes = size(xx,1);

% GDof: global number of degrees of freedom
GDof = 2*numberNodes;

% computation of the system stiffness matrix
[stiffness,force,mass] = ...
    formStiffnessMassTimoshenkoBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,C,0,rho,I,thickness);

% boundary conditions (simply-supported at both bords)
% fixedNodeW = [1 ; numberNodes];
% fixedNodeTX = [];
% boundary conditions (clamped at both bords)
fixedNodeW = [1 ; numberNodes];
fixedNodeTX = fixedNodeW;
% boundary conditions (cantilever)
% fixedNodeW = [1];
% fixedNodeTX = [1];
prescribedDof = [fixedNodeW; fixedNodeTX+numberNodes];

% Time transient simulation
timeStep = 0.005;
totalTime = 0.5;
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
ws = 1:numberNodes;
displacementsTime(ws,1) = sin(pi*xx) - pi*xx.*(1-xx);
displacementsTime(ws+numberNodes,1) = -pi*cos(pi*xx) + pi*(1-2*xx);
accelerationsTime(:,1) = mass\(force ...
    - stiffness*displacementsTime(:,1));

for i = 2:length(time)
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

% central point vs time
figure;
plot(time,displacementsTime(round(numberNodes/2),:),...
    '.-','linewidth',2,'markersize',16)
xlabel('time'); ylabel('central point motion')
ylim([-0.24 0.24])
set(gca,'linewidth',2,'fontsize',14)
grid on; box on;
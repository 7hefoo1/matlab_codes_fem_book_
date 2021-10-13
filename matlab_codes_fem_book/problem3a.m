% ................................................................
% MATLAB codes for Finite Element Analysis
% problem3a.m
% ref: D. Logan, A first course in the finite element method,
% third Edition, page 121, exercise P3-10
% with isoparametric formulation
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E = 70000; A = 200; EA = E*A; k = 2000; 

% generation of coordinates and connectivities
numberElements = 3;  
numberNodes = 4;
elementNodes = [1 2; 2 3; 3 4];
nodeCoordinates = [0 2000 4000 4000];
xx = nodeCoordinates;

% for structure:
%   displacements: displacement vector
%   force : force vector
%   stiffness: stiffness matrix
displacements = zeros(numberNodes,1);
force = zeros(numberNodes,1);
stiffness = zeros(numberNodes,numberNodes); 

% applied load at node 2
force(2) = 8000.0;

% computation of the system stiffness matrix
ea = zeros(1,numberElements);
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    elementDof = elementNodes(e,:);
    
    if e < 3 % bar elements
        nn = length(elementDof);
        length_element = nodeCoordinates(elementDof(2)) ...
                         -nodeCoordinates(elementDof(1));
        detJacobian = length_element/2;
        invJacobian = 1/detJacobian;
        % central Gauss point (xi=0, weight W=2)
        [shape,naturalDerivatives] = shapeFunctionL2(0.0);
        Xderivatives = naturalDerivatives*invJacobian;
        
        % B matrix
        B = zeros(1,nn); B(1:nn) = Xderivatives(:);
        ea(e) = E*A;
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) ...
            + B'*B*2*detJacobian*ea(e);
    else % spring element
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) + k*[1 -1;-1 1];
    end
end

% boundary conditions and solution
prescribedDof = [1;4];

% solution
displacements = solution(numberNodes,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    numberNodes,prescribedDof)
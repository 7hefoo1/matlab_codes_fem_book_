% ................................................................
% MATLAB codes for Finite Element Analysis
% problem3vib.m
% ref: J.N. Reddy, An introduction to the Finite Element Method,
% third Edition, page 86, example 2.5.4
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
% rho: density
E = 70000; A = 200; EA = E*A; k = EA/4000; rho = 1000;

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
%   mass: mass matrix
displacements = zeros(numberNodes,1);
force = zeros(numberNodes,1);
stiffness = zeros(numberNodes,numberNodes); 
mass = zeros(numberNodes,numberNodes); 

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
        
        % stiffness matrix. Central Gauss point (xi=0, weight W=2)
        [shape,naturalDerivatives] = shapeFunctionL2(0.0);
        Xderivatives = naturalDerivatives*invJacobian;
        
        % B matrix
        B = zeros(1,nn); B(1:nn) = Xderivatives(:);
        ea(e) = E*A;
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) ...
            + B'*B*2*detJacobian*ea(e);
        
        % mass matrix
        % exact integration -------
%         mass(elementDof,elementDof) = ...
%             mass(elementDof,elementDof) ...
%             + [2 1;1 2]*detJacobian*rho*A/3;
        % lumped mass matrix -------
%         mass(elementDof,elementDof) = ...
%             mass(elementDof,elementDof) ...
%             + [1 0;0 1]*detJacobian*rho*A;
        % Gauss quadrature calculation -------
        % two-points integration (coincides with exact)
%         gaussLocations = [0.577350269189626; -0.577350269189626];
%         gaussWeights = ones(2,1);
        % one-point integration (reduced integration)
        gaussLocations = 0.0;
        gaussWeights = 2;
        for q = 1:size(gaussWeights,1)
            [shape,~] = shapeFunctionL2(gaussLocations(q));
            
            mass(elementDof,elementDof) = ...
                mass(elementDof,elementDof) ...
                + shape*shape'*gaussWeights(q)*detJacobian*rho*A;
        end
    else % spring element
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) + k*[1 -1;-1 1];
    end
end

% boundary conditions and solution
prescribedDof = [1;4];

% free vibration problem
[modes,eigenvalues] = eigenvalue(numberNodes,prescribedDof,...
                                 stiffness,mass,0);

omega = sqrt(eigenvalues)*sqrt(rho/E)*4000
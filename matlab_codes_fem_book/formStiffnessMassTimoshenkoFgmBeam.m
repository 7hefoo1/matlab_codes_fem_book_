function [stiffness,force,mass] = ...
    formStiffnessMassTimoshenkoFgmBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,C,P,I,thickness)

% computation of stiffness, mass matrices and force 
% vector for Timoshenko beam element
stiffness = zeros(GDof);
mass = zeros(GDof);
force = zeros(GDof,1);

% 2x2 Gauss quadrature
gaussLocations = [0.577350269189626;-0.577350269189626];
gaussWeights = ones(2,1);

% bending contribution for matrices
for e = 1:numberElements
    indice = elementNodes(e,:);
    elementDof = [indice indice+numberNodes indice+2*numberNodes];
    indiceMass = indice+numberNodes;
    
    ndof = length(indice);
    length_element = xx(indice(2))-xx(indice(1));
    detJacobian = length_element/2; invJacobian=1/detJacobian;
    for q = 1:size(gaussWeights,1)
        pt = gaussLocations(q,:);
        [shape,naturalDerivatives] = shapeFunctionL2(pt(1));
        Xderivatives = naturalDerivatives*invJacobian;
        % B matrix
        B = zeros(3,3*ndof);
        B(1,1:ndof) = Xderivatives(:)';
        B(2,2*ndof+1:3*ndof) = Xderivatives(:)';
        
        % stiffness matrix
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) + ...
            B'*C*B*gaussWeights(q)*detJacobian;
        
        % force vector
        force(indiceMass) = force(indiceMass) + ...
            shape*P*detJacobian*gaussWeights(q);
        
        % B matrix
        B = zeros(3,3*ndof);
        B(1,1:ndof) = shape;
        B(2,1+ndof:2*ndof) = shape;
        B(3,1+2*ndof:3*ndof) = shape;
        
        % mass matrix
        mass(elementDof,elementDof) = ...
            mass(elementDof,elementDof) + ...
            B'*I*B*gaussWeights(q)*detJacobian;
    end
end

% shear contribution for the matrices
gaussLocations = 0.;
gaussWeights = 2.;

for e = 1:numberElements
    indice = elementNodes(e,:);
    elementDof = [ indice indice+numberNodes indice+2*numberNodes];
    ndof = length(indice);
    length_element = xx(indice(2))-xx(indice(1));
    detJ0 = length_element/2; invJ0 = 1/detJ0;
    for q = 1:size(gaussWeights,1)
        pt = gaussLocations(q,:);
        [shape,naturalDerivatives] = shapeFunctionL2(pt(1));
        Xderivatives = naturalDerivatives*invJacobian;
        % B
        B = zeros(3,3*ndof);
        B(3,ndof+1:2*ndof) = Xderivatives(:)';
        B(3,2*ndof+1:3*ndof) = shape;
        
        % stiffness matrix
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) + ...
            B'*C*B*gaussWeights(q)*detJacobian;
    end
end

end
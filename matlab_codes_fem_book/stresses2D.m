function [stress,strain] = stresses2D(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates, ...
    displacements,C,elemType,quadType)

% quadrature according to quadType
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);

% stresses at nodes
stress = zeros(numberElements,size(gaussLocations,1),3);
% stressPoints = [-1 -1;1 -1;1 1;-1 1];

for e = 1:numberElements
    indice = elementNodes(e,:);
    elementDof = [ indice indice+numberNodes ];
    nn = length(indice);
    for q = 1:size(gaussWeights,1)
        pt = gaussLocations(q,:);
        wt = gaussWeights(q);
        xi = pt(1);
        eta = pt(2);
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

        %  B matrix
        B = zeros(3,2*nn);
        B(1,1:nn)       = XYderivatives(:,1)';
        B(2,nn+1:2*nn)  = XYderivatives(:,2)';
        B(3,1:nn)       = XYderivatives(:,2)';
        B(3,nn+1:2*nn)  = XYderivatives(:,1)';
        
        % element deformation
        strain = B*displacements(elementDof);
        stress(e,q,:) = C*strain;
    end
end

end  % end function stresses2D
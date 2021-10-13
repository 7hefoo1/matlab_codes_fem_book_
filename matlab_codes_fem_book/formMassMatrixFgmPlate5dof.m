function [M] = ...
    formMassMatrixFgmPlate5dof(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,I,...
    elemType,quadType)

% computation of mass matrix for Mindlin plate element

M = zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);

% cycle for element
for e=1:numberElements
    % indice: nodal connectivities for each element
    % elementDof: element degrees of freedom
    indice=elementNodes(e,:);
    elementDof = [indice indice+numberNodes indice+2*numberNodes ...
        indice+3*numberNodes indice+4*numberNodes];
    ndof=length(indice);
    
    % cycle for Gauss point
    for q=1:size(gaussWeights,1)
        GaussPoint=gaussLocations(q,:);
        xi=GaussPoint(1);
        eta=GaussPoint(2);
        
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % [N] matrix
        N = zeros(5,5*ndof);
        N(1,1:ndof)          = shapeFunction;
        N(2,ndof+1:2*ndof)   = shapeFunction;
        N(3,2*ndof+1:3*ndof) = shapeFunction;
        N(4,3*ndof+1:4*ndof) = shapeFunction;
        N(5,4*ndof+1:5*ndof) = shapeFunction;
        
        % mass matrix
        M(elementDof,elementDof) = M(elementDof,elementDof) + ...
            N'*I*N*gaussWeights(q)*det(Jacob);
        
        
    end  % end Gauss point loop
end % end element loop

end
function [force] = ...
    formForceVectorK(GDof,numberElements,elementNodes,...
    numberNodes,nodeCoordinates,P,quadType,dof_per_node)

% computation of force vector for Kirchhoff plate element

% force: force vector
force = zeros(GDof,1);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);

% cycle for element
for e = 1:numberElements
    % indice : nodal connectivities for each element
    indice = elementNodes(e,:);
    if dof_per_node == 3
        elementDof = [indice indice+numberNodes ...
            indice+2*numberNodes];
    else % 4 dof
        elementDof = [indice indice+numberNodes ...
            indice+2*numberNodes indice+3*numberNodes];
    end
    ndof = length(elementDof);
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(q,:);
        GaussWeight = gaussWeights(q);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % part related to the mapping
        % shape functions and derivatives
        [~,natDerQ4] = shapeFunctionKQ4(xi,eta);
        if dof_per_node == 3
            [shapeFunction,~] = shapeFunctionK12(xi,eta);
        else % 4 dof
            [shapeFunction,~] = shapeFunctionK16(xi,eta);
        end
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,~,~] = JacobianK(nodeCoordinates(indice,:),natDerQ4);
        
        % force vector
        force(elementDof) = force(elementDof) + ...
            shapeFunction*P*det(Jacob)*GaussWeight;
    end  % end Gauss point loop
    
end % end element loop

end

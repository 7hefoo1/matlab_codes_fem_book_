function [JacobianMatrix,invJacobian,XYDerivatives] = ...
    JacobianK(nodeCoordinates,naturalDerivatives)

% JacobianMatrix: Jacobian matrix
% invJacobian: inverse of Jacobian Matrix
% XYDerivatives: derivatives w.r.t. x and y
% naturalDerivatives: derivatives w.r.t. xi and eta
% nodeCoordinates: nodal coordinates at element level

JacobianMatrix = nodeCoordinates'*naturalDerivatives(:,1:2);
invJacobian = inv(JacobianMatrix);

XYDerivatives = nodeCoordinates'*naturalDerivatives;

end % end function Jacobian

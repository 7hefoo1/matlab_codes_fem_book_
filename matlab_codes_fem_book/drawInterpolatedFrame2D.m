% reordering displacements and rotations
U = displacements(1:numberNodes);
W = displacements(numberNodes+1:2*numberNodes);
R = displacements(2*numberNodes+1:3*numberNodes);


% graphical representation with interpolation for each element
interpNodes = 50;
hold on; box on;
for n = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(n,:);
    elementDof = [indice indice+numberNodes indice+2*numberNodes];
    nn = length(indice);
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    cosa = xa/length_element;
    sena = ya/length_element;
    ll = length_element;
    
    L = [cosa*eye(2) sena*eye(2);
        -sena*eye(2) cosa*eye(2)];
    
    % displacements in local coordinates
    D = L*[U(indice(1));U(indice(2));W(indice(1));W(indice(2))];
    
    % interpolated local coordinates
    nodeInterp(:,1) = linspace(nodeCoordinates(indice(1),1), ...
        nodeCoordinates(indice(2),1),interpNodes);
    nodeInterp(:,2) = linspace(nodeCoordinates(indice(1),2), ...
        nodeCoordinates(indice(2),2), interpNodes);
    % local coordinate
    if cosa == 0 % beam is vertical
        xp = abs(nodeInterp(:,2)-nodeInterp(1,2));
    else
        xp = abs(nodeInterp(:,1)-nodeInterp(1,1))/cosa;
    end
    % Lagrange shape function (stretching)
    psi1 = 1 - xp/ll;
    psi2 = xp/ll;
    % Hermite shape function (bending)
    xi = 2*xp/ll-1;
    phi1 = 0.25*(2 - 3*xi + xi.^3);
    phi2 = ll/8*(1 - xi - xi.^2 + xi.^3);
    phi3 = 0.25*(2 + 3*xi - xi.^3);
    phi4 = ll/8*(-1 - xi + xi.^2 + xi.^3);
        
    u = psi1*D(1) + psi2*D(2);
    w = phi1*D(3) + phi2*R(indice(1)) + ...
        phi3*D(4) + phi4*R(indice(2));
    L = [cosa*eye(interpNodes) -sena*eye(interpNodes);
        sena*eye(interpNodes) cosa*eye(interpNodes)];
    DInterp = L*[u;w];
    elementNodesInterp = [[1:interpNodes-1]' [2:interpNodes]'];
    drawingMesh(nodeInterp+scaleFact*[DInterp(1:interpNodes) DInterp(interpNodes+1:2*interpNodes)],...
        elementNodesInterp, 'L2','k-');
end
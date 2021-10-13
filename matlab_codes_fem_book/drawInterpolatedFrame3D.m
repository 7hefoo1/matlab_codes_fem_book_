% reordering displacements and rotations
ux = U(1:6:6*numberNodes);
uy = U(2:6:6*numberNodes);
uz = U(3:6:6*numberNodes);
rx = U(4:6:6*numberNodes);
ry = U(5:6:6*numberNodes);
rz = U(6:6:6*numberNodes);

% graphical representation with interpolation for each element
interpNodes = 50;
hold on; box on;
for n = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(n,:);
    elementDof = [6*indice(1)-5 6*indice(1)-4 6*indice(1)-3 ...
        6*indice(1)-2 6*indice(1)-1 6*indice(1)...
        6*indice(2)-5 6*indice(2)-4 6*indice(2)-3 ...
        6*indice(2)-2 6*indice(2)-1 6*indice(2)] ;
    x1 = nodeCoordinates(indice(1),1);
    y1 = nodeCoordinates(indice(1),2);
    z1 = nodeCoordinates(indice(1),3);
    x2 = nodeCoordinates(indice(2),1);
    y2 = nodeCoordinates(indice(2),2);
    z2 = nodeCoordinates(indice(2),3);
    
    L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + ...
        (z2-z1)*(z2-z1));
    
    if x1 == x2 && y1 == y2
        if z2 > z1
            Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
        else
            Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
        end
    else
        CXx = (x2-x1)/L;
        CYx = (y2-y1)/L;
        CZx = (z2-z1)/L;
        D = sqrt(CXx*CXx + CYx*CYx);
        CXy = -CYx/D;
        CYy = CXx/D;
        CZy = 0;
        CXz = -CXx*CZx/D;
        CYz = -CYx*CZx/D;
        CZz = D;
        Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
        
    end
    R = [Lambda zeros(3,9); zeros(3) Lambda zeros(3,6);
        zeros(3,6) Lambda zeros(3);zeros(3,9) Lambda];
    
    % displacements in local coordinates
    d = R*[ux(indice(1));uy(indice(1));uz(indice(1));
        rx(indice(1));ry(indice(1));rz(indice(1));
        ux(indice(2));uy(indice(2));uz(indice(2));
        rx(indice(2));ry(indice(2));rz(indice(2))];
    
    % interpolated local coordinates
    nodeInterp(:,1) = linspace(nodeCoordinates(indice(1),1), ...
        nodeCoordinates(indice(2),1),interpNodes);
    nodeInterp(:,2) = linspace(nodeCoordinates(indice(1),2), ...
        nodeCoordinates(indice(2),2), interpNodes);
    nodeInterp(:,3) = linspace(nodeCoordinates(indice(1),3), ...
        nodeCoordinates(indice(2),3), interpNodes);
    % local coordinate
    if x1 == x2 && y1 == y2 % beam along z
        xp = abs(nodeInterp(:,3)-nodeInterp(1,3));
    elseif CXx == 0 % beam along y
        xp = abs(nodeInterp(:,2)-nodeInterp(1,2));
    else
        xp = abs(nodeInterp(:,1)-nodeInterp(1,1))/abs(CXx);
    end
    % Lagrange shape function (stretching)
    psi1 = 1 - xp/L;
    psi2 = xp/L;
    % Hermite shape function (bending)
    xi = 2*xp/L-1;
    phi1 = 0.25*(2 - 3*xi + xi.^3);
    phi2 = L/8*(1 - xi - xi.^2 + xi.^3);
    phi3 = 0.25*(2 + 3*xi - xi.^3);
    phi4 = L/8*(-1 - xi + xi.^2 + xi.^3);
        
    u = psi1*d(1) + psi2*d(6+1);
    v = phi1*d(2) + phi2*d(6) + ...
        phi3*d(6+2) + phi4*d(6+6);
    w = phi1*d(3) + phi2*d(5) + ...
        phi3*d(6+3) + phi4*d(6+5);
    for i = 1:interpNodes
        DInterp(1+(i-1)*3:i*3,1) = Lambda\[u(i);v(i);w(i)];
    end
    
    elementNodesInterp = [[1:interpNodes-1]' [2:interpNodes]'];
    drawingMesh(nodeInterp+scaleFact*[DInterp(1:3:3*interpNodes) DInterp(2:3:3*interpNodes) DInterp(3:3:3*interpNodes)],...
        elementNodesInterp, 'L2','k-');
end
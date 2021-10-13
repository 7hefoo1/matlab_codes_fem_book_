interpNodes = 50;
hold on; box on;
for n = 1:numberElements
    nodeA = elementNodes(n,1); nodeB = elementNodes(n,2);
    XX = linspace(nodeCoordinates(nodeA),nodeCoordinates(nodeB),...
        interpNodes);
    ll = nodeCoordinates(nodeB)-nodeCoordinates(nodeA);
    % dimensionless shape functions in Cartesian coordinates
    xi = (XX - nodeCoordinates(nodeA))*2/ll - 1;
    % Hermite shape function
    phi1 = 0.25*(2 - 3*xi + xi.^3);
    phi2 = ll/8*(1 - xi - xi.^2 + xi.^3);
    phi3 = 0.25*(2 + 3*xi - xi.^3);
    phi4 = ll/8*(-1 - xi + xi.^2 + xi.^3);
    
    w = phi1*W(nodeA) + phi2*R(nodeA) + phi3*W(nodeB) + phi4*R(nodeB);
    plot(XX,w,'-k','linewidth',1.5)
end
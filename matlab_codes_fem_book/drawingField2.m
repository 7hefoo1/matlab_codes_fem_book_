function drawingField2(nodeCoordinates,elementNodes,...
    scaleFactor,UX,UY,elementType,field)

switch elementType
    case 'Q4'
        for i = 1:size(elementNodes,1)
            XX = nodeCoordinates(elementNodes(i,1:4),1) + scaleFactor*UX(elementNodes(i,1:4));
            YY = nodeCoordinates(elementNodes(i,1:4),2) + scaleFactor*UY(elementNodes(i,1:4));
            patch(XX,YY,field(i,:))
        end
    case {'Q8', 'Q9'}
        for i = 1:size(elementNodes,1)
            XX = nodeCoordinates(elementNodes(i,1:4),1) + scaleFactor*UX(elementNodes(i,:));
            YY = nodeCoordinates(elementNodes(i,1:4),2) + scaleFactor*UY(elementNodes(i,:));
            patch(XX,YY,field(i,:))
        end
end


end



function drawingMesh(nodeCoordinates,elementNodes,elementType,lineType)
%drawingMesh mesh representation for 1D and 2D FEs
% Input data required
%
% nodeCoordinates : Cartesian coordinates of the nodes
% elementNodes    : Element connectivity
% elementType     : Type of element (e.g. linear, quadratic, etc.)
% lineType        : Type of line for plotting purposes
%

%%

% resize the vector to suit 2d and 3d problems (beams)
if size(nodeCoordinates,2) == 2
    nodeCoordinates(end,3) = 0;
end

switch elementType
    case 'L2' % linear element in 2D/3D with 2 nodes
        for k = 1:size(elementNodes,1)
            plot3(nodeCoordinates(elementNodes(k,:),1),nodeCoordinates(elementNodes(k,:),2),nodeCoordinates(elementNodes(k,:),3),lineType,'MarkerSize',12,'Linewidth',1.5)
        end
        
    case 'L3' % linear element in 2D/3D with 3 nodes
        for k = 1:size(elementNodes,1)
            plot3(nodeCoordinates(elementNodes(k,:),1),nodeCoordinates(elementNodes(k,:),2),nodeCoordinates(elementNodes(k,:),3),lineType,'MarkerSize',12,'Linewidth',1.5)
        end
        
    case 'Q4'
        for k = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(k,:),1),nodeCoordinates(elementNodes(k,:),2),'w','FaceColor','none','LineStyle',lineType,'EdgeColor','k')
        end

    case {'Q8','Q9'}
        for k = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(k,1:4),1),nodeCoordinates(elementNodes(k,1:4),2),'w','FaceColor','none','LineStyle',lineType,'EdgeColor','k')
        end
    otherwise
        disp('Element type not available')
end

end


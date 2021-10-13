function [node, element] = ...
    rectangularMesh(Lx, Ly, nelemX, nelemY, ...
    elementType)

% element size
deltaX = Lx/nelemX;
deltaY = Ly/nelemY;

switch elementType
    case 'Q4' % Q4 element
        
        nodesX = nelemX+1;
        nodesY = nelemY+1;
        
        % nodal coordinates
        node = [];
        for j = 1:nodesY
            for i = 1:nodesX
                x = (i-1)*deltaX; y = (j-1)*deltaY;
                node = [node; x y];
            end
        end
        
        % connectivity
        element = [];
        for j = 1:nelemY
            for i = 1:nelemX
                i1 = i+(j-1)*nodesX;
                i2 = i1+1;
                i3 = i2+nodesX;
                i4 = i1+nodesX;
                element = [element; i1 i2 i3 i4];
            end
        end
        
    case 'Q8' % Q8 element
        nodesX = nelemX*2+1;
        nodesY = nelemY*2+1;
        
        % nodal coordinates
        node = [];
        for j = 1:nodesY
            y = (j-1)*deltaY/2;
            if mod(j,2) % odd index
                nodesX = nelemX*2+1;
                for i = 1:nodesX
                    x = (i-1)*deltaX/2;
                    node = [node; x y];
                end
            else % even index
                nodesX = nelemX+1;
                for i = 1:nodesX
                    x = (i-1)*deltaX;
                    node = [node; x y];
                end
            end
        end
        
        % connectivity
        element = [];
        for j = 1:nelemY
            for i = 1:nelemX
                i1 = 1 + 2*(i-1) + (j-1)*round(1.5*nodesX);
                i2 = i1 + 2;
                i3 = i2 + round(1.5*nodesX);
                i4 = i1 + round(1.5*nodesX);
                i5 = i1 + 1;
                i6 = i1 + nodesX+2-i;
                i7 = i4 + 1;
                i8 = i6 - 1;
                element = [element; i1 i2 i3 i4 i5 i6 i7 i8];
            end
        end

    case 'Q9' % Q9 element        
        nodesX = nelemX*2+1;
        nodesY = nelemY*2+1;
        
        % nodal coordinates
        node = [];
        for j = 1:nodesY
            for i = 1:nodesX
                x = (i-1)*deltaX/2;
                y = (j-1)*deltaY/2;
                node = [node; x y];
            end
        end
        
        % connectivity
        element = [];
        for j = 1:nelemY
            for i = 1:nelemX
                i1 = 1+2*(i-1)+(j-1)*2*nodesX;
                i2 = i1 + 2;
                i3 = i2 + 2*nodesX;
                i4 = i1 + 2*nodesX;
                i5 = i1 + 1;
                i6 = i2 + nodesX;
                i7 = i4 + 1;
                i8 = i1 + nodesX;
                i9 = i5 + nodesX;
                element = [element; i1 i2 i3 i4 i5 i6 i7 i8 i9];
            end
        end
end

end
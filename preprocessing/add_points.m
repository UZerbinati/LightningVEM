function [domainMesh] = add_points(domainMesh, n)

n1   = n + 1;
init = domainMesh.nvertex;

check = sparse(domainMesh.nvertex, domainMesh.nvertex);

for i = 1:domainMesh.npolygon

    number  = numel(domainMesh.connect{i});
    newconn = zeros(number * n1, 1);

    for j = 1 : number-1
        
        newconn((j-1)*n1 + 1) = domainMesh.connect{i}(j);  

        p1 = domainMesh.connect{i}(j);
        p2 = domainMesh.connect{i}(j+1);

        a = domainMesh.coords(domainMesh.connect{i}(j),:);
        b = domainMesh.coords(domainMesh.connect{i}(j+1),:);
        
        if (check(p1,p2) == 0)

            check(p1,p2) = init;
            check(p2,p1) = init;

            xnodes = linspace(a(1),b(1),n+2);
            ynodes = linspace(a(2),b(2),n+2);
    
            domainMesh.coords(init+1:init+n,1) = xnodes(2:end-1);
            domainMesh.coords(init+1:init+n,2) = ynodes(2:end-1);
    
            for k = 1:n
                init = init + 1;
                newconn((j-1)*n1 + 1 + k) = init;
    
                if(domainMesh.coords(init,1) == 0)
                    domainMesh.boundary_nodes.left = [domainMesh.boundary_nodes.left; init];
                    domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
                end
                if(domainMesh.coords(init,1) == 1)
                    domainMesh.boundary_nodes.right = [domainMesh.boundary_nodes.right; init];
                    domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
                end
                if(domainMesh.coords(init,2) == 0)
                    domainMesh.boundary_nodes.bottom = [domainMesh.boundary_nodes.bottom; init];
                    domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
                end
                if(domainMesh.coords(init,2) == 1)
                    domainMesh.boundary_nodes.top = [domainMesh.boundary_nodes.top; init];
                    domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
                end

            end

        else

            for k = 1:n

                newconn((j-1)*n1 + 1 + k) = check(p1,p2) + (n+1-k);

            end
        
        end

    end

    newconn((number-1)*n1 + 1) = domainMesh.connect{i}(number);

    p1 = domainMesh.connect{i}(end);
    p2 = domainMesh.connect{i}(1);

    a = domainMesh.coords(domainMesh.connect{i}(end),:);
    b = domainMesh.coords(domainMesh.connect{i}(1),:);

    if (check(p1,p2) == 0)

        check(p1,p2) = init;
        check(p2,p1) = init;
    
        xnodes = linspace(a(1),b(1),n+2);
        ynodes = linspace(a(2),b(2),n+2);
       
        domainMesh.coords(init+1:init+n,1) = xnodes(2:end-1);
        domainMesh.coords(init+1:init+n,2) = ynodes(2:end-1);
    
        for k = 1:n
            init = init + 1;
            newconn((number-1)*n1 + 1 + k) = init;
    
            if(domainMesh.coords(init,1) == 0)
                domainMesh.boundary_nodes.left = [domainMesh.boundary_nodes.left; init];
                domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
            end
            if(domainMesh.coords(init,1) == 1)
                domainMesh.boundary_nodes.right = [domainMesh.boundary_nodes.right; init];
                domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
            end
            if(domainMesh.coords(init,2) == 0)
                domainMesh.boundary_nodes.bottom = [domainMesh.boundary_nodes.bottom; init];
                domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
            end
            if(domainMesh.coords(init,2) == 1)
                domainMesh.boundary_nodes.top = [domainMesh.boundary_nodes.top; init];
                domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.all; init];
            end
        end

    else

        for k = 1:n

            newconn((number-1)*n1 + 1 + k) = check(p1,p2) + (n+1-k);

        end

    end
    
    domainMesh.connect{i} = newconn;
end

domainMesh.nvertex = init;

end
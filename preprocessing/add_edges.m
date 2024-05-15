function [domainMesh] = add_edges(domainMesh,k) 

% add_edges: this function adds the information relative to the edges to the struct domainMesh
%
% Input parameters:
% domainMesh: the struct relative to the mesh of which we want the information;
%          k: degree of the polynomials that will be used.
%
% Output parameters:
% domainMesh: the struct now contains also the following fields:
%        adj: each row contains the two elements adjacent to the edge
%      edges: a sparse upper triangular matrix in  which we store in  (i,j) the index of the 
%             edge joinin i and j (if it exists)
%     length: the lenghts of the edges
%        map: a matrix in which we store in the i-th row the vertices joined by the i-th edge
%     nedges: the number of edges

domainMesh.edges     = sparse(domainMesh.nvertex, domainMesh.nvertex);                               %Initialization
domainMesh.map       = zeros(4 * domainMesh.nvertex, 2);                                             %Number of edges is usually lower than 4 * nvertex
domainMesh.adj       = zeros(4 * domainMesh.nvertex, 2);
domainMesh.length    = zeros(4 * domainMesh.nvertex, 1);
domainMesh.intcoords = zeros(4 * (k - 1) * domainMesh.nvertex, 2);
counter = 1;                                                                                         %Edges counter 

lobatto = (lobatto_points_1d(k+1) + 1) / 2;

for i=1:domainMesh.npolygon                                                                          %For each polygon in the mesh...
    
    for j=1:(numel(domainMesh.connect{i,1})-1)                              
        
        a = domainMesh.connect{i,1}(j);                                                              %j-th   vertex
        b = domainMesh.connect{i,1}(j+1);                                                            %j+1-th vertex

        if (a > b)                                                                                   %We want the matrix to be upper triangular. If a is bigger
                                                                                                     %than b we switch the two indices            
            temp = a;
            a = b;
            b = temp;
            
        end

        xa = domainMesh.coords(a,1);
        xb = domainMesh.coords(b,1);
        ya = domainMesh.coords(a,2);
        yb = domainMesh.coords(b,2);
        
        
        if (domainMesh.edges(a,b) == 0)                                                              %Check if the edge doesn't exist yet
            
            domainMesh.edges(a,b)      = counter;                                                    %Add the edge
            domainMesh.map(counter,:)  = [a b];                                                      %Store the vertices joined by the edge
            domainMesh.adj(counter,1)  = i;                                                          %Add the first element that has the edge
            domainMesh.length(counter) = sqrt((xa-xb)^2 + (ya-yb)^2);                                %Compute the length
            
            if (k > 1)
                
                h     = domainMesh.length(counter) / k;
                index = ((counter-1)*(k-1) + 1 : counter*(k-1));

                domainMesh.intcoords(index,1) = xa + lobatto(2:end-1)*(xb-xa);
                domainMesh.intcoords(index,2) = ya + lobatto(2:end-1)*(yb-ya);

            end

            counter = counter + 1;                                                                   %Update the counter

        else 
         
            domainMesh.adj(domainMesh.edges(a,b),2) = i;                                             %Add the second element that has the edge
         
        end
    end
    
    a = domainMesh.connect{i,1}(end);                                                                %Analogous for the last and the first vertex
    b = domainMesh.connect{i,1}(1);

    if (b<a)
            
         temp = a;
         a = b;
         b = temp;
            
    end

    xa = domainMesh.coords(a,1);
    xb = domainMesh.coords(b,1);
    ya = domainMesh.coords(a,2);
    yb = domainMesh.coords(b,2);
        
    
     if (domainMesh.edges(a,b) == 0)
            
         domainMesh.edges(a,b)      = counter;
         domainMesh.map(counter,:)  = [a b];
         domainMesh.adj(counter,1)  = i;
         domainMesh.length(counter) = sqrt((xa-xb)^2 + (ya-yb)^2);

         if (k > 1)
                
            index = ((counter-1)*(k-1) + 1 : counter*(k-1));

            domainMesh.intcoords(index,1) = xa + lobatto(2:end-1)*(xb-xa);
            domainMesh.intcoords(index,2) = ya + lobatto(2:end-1)*(yb-ya);

        end

         counter = counter + 1;
     
     else 
         
         domainMesh.adj(domainMesh.edges(a,b),2) = i;
         
     end
    
    
end

domainMesh.map       = domainMesh.map(1:counter-1,:);                                                %Remove the unusued rows
domainMesh.adj       = domainMesh.adj(1:counter-1,:);
domainMesh.length    = domainMesh.length(1:counter-1);
domainMesh.nedges    = counter - 1;
domainMesh.intcoords = domainMesh.intcoords(1:(counter-1)*(k-1),:);

end
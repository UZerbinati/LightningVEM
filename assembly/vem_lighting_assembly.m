function [K_global,f_global,domainMesh] = vem_lighting_assembly(domainMesh, matProps, polynomial, f, k)

%% MEMORY ALLOCATION
siz = domainMesh.nvertex + (polynomial.k - 1) * domainMesh.nedges ...                                %Size of the linear system
    + polynomial.int * domainMesh.npolygon;      

K_global = zeros(siz,siz);                                                                           %Memory allocation
f_global = zeros(siz, 1);                                                                          

domainMesh.polygon  = {1:domainMesh.npolygon};                                                       %Allocate memory for polygon information
domainMesh.diameter = 0;                                                                             %Initialize diameter of the mesh

connect  = domainMesh.connect;
coords   = parallel.pool.Constant(domainMesh.coords);
polygons = cell(domainMesh.npolygon,1);
K_loc    = cell(domainMesh.npolygon,1);
f_loc    = cell(domainMesh.npolygon,1);

%% CONSTRUCT THE LINEAR SYTEM
parfor i = 1:domainMesh.npolygon
   
    local_vertex = connect{i,1};                                                                     %Get the vertex indexes of the i-th polygon

    vertex = coords.Value(local_vertex,:);

    [K_loc{i}, f_loc{i}, polygons{i}] = vem_lighting_element(vertex, matProps, f, k);             %Assembly local stiffness matrix and load term
    
end

for i = 1:domainMesh.npolygon
    
    local_vertex = connect{i,1};
    %local_edges  = get_edges(domainMesh.edges, connect{i,1});                                       %Get the edge indexes of the i-th polygon
    %local_dofs   = get_local_dofs(domainMesh, polynomial, local_vertex,local_edges);                %CHECK!!

    polygons{i}.local_dofs = local_vertex;                                                           %Store the polygon information 
    
    if (polygons{i}.diameter > domainMesh.diameter)                                                  %Update the diameter of the mesh

        domainMesh.diameter = polygons{i}.diameter;

    end

    K_global(local_vertex, local_vertex) = K_global(local_vertex, local_vertex) + K_loc{i};
    f_global(local_vertex)               = f_global(local_vertex)               + f_loc{i};

end

domainMesh.polygon = polygons;
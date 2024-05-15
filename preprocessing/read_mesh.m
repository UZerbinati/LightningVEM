function domainMesh = read_mesh(meshFile)

%% READ DOMAIN TYPE 
mesh_file   = fopen(meshFile);                                              %Open mesh file
line        = fgets(mesh_file);                                             %Read commented line  
domain_type = sscanf(fgets(mesh_file),'%s');                                %Read domain type
fclose(mesh_file);    
  
if     strcmp(domain_type,'RectangularDomain')                              %Read the mesh
    domainMesh = read_mesh_rectangular_domain_poisson2d(meshFile);
elseif strcmp(domain_type,'WrenchDomain')
    domainMesh = read_mesh_wrench_domain_poisson2d(meshFile);      
elseif strcmp(domain_type,'PlateWithHoleDomain')
    domainMesh = read_mesh_plate_with_hole_domain_poisson2d(meshFile);       
else
    throw_error('In read_mesh.m: domain_type\n');      
end    

function domainMesh = read_mesh_rectangular_domain_linelast2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the bottom boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.bottom = sscanf(line,'%d');
  % read indices of nodes that are located on the top boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.top = sscanf(line,'%d');
  % read indices of nodes that are located on the left boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.left = sscanf(line,'%d');    
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.right = sscanf(line,'%d');  
  % all boundary nodes
  domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.bottom;domainMesh.boundary_nodes.top;...
                             domainMesh.boundary_nodes.left;domainMesh.boundary_nodes.right];
  domainMesh.boundary_nodes.all=unique(domainMesh.boundary_nodes.all); % delete repeated nodes from boundary_nodes.all
  % read xmax, xmin, ymax, ymin for the rectangular domain
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs
  % all
  domainMesh.boundary_dofs.all=zeros(2*length(domainMesh.boundary_nodes.all),1);      
  range=1:length(domainMesh.boundary_nodes.all);
  domainMesh.boundary_dofs.all(2*range-1)=2*domainMesh.boundary_nodes.all-1;
  domainMesh.boundary_dofs.all(2*range)=2*domainMesh.boundary_nodes.all;    
  % bottom
  domainMesh.boundary_dofs.bottom=zeros(2*length(domainMesh.boundary_nodes.bottom),1);      
  range=1:length(domainMesh.boundary_nodes.bottom);
  domainMesh.boundary_dofs.bottom(2*range-1)=2*domainMesh.boundary_nodes.bottom-1;
  domainMesh.boundary_dofs.bottom(2*range)=2*domainMesh.boundary_nodes.bottom;    
  % top
  domainMesh.boundary_dofs.top=zeros(2*length(domainMesh.boundary_nodes.top),1);      
  range=1:length(domainMesh.boundary_nodes.top);
  domainMesh.boundary_dofs.top(2*range-1)=2*domainMesh.boundary_nodes.top-1;
  domainMesh.boundary_dofs.top(2*range)=2*domainMesh.boundary_nodes.top;   
  % left
  domainMesh.boundary_dofs.left=zeros(2*length(domainMesh.boundary_nodes.left),1);      
  range=1:length(domainMesh.boundary_nodes.left);
  domainMesh.boundary_dofs.left(2*range-1)=2*domainMesh.boundary_nodes.left-1;
  domainMesh.boundary_dofs.left(2*range)=2*domainMesh.boundary_nodes.left;   
  % right
  domainMesh.boundary_dofs.right=zeros(2*length(domainMesh.boundary_nodes.right),1);      
  range=1:length(domainMesh.boundary_nodes.right);
  domainMesh.boundary_dofs.right(2*range-1)=2*domainMesh.boundary_nodes.right-1;
  domainMesh.boundary_dofs.right(2*range)=2*domainMesh.boundary_nodes.right;   
  domainMesh.nvertex  = numel(domainMesh.coords) / 2;                                                  %Number of vertices
  domainMesh.npolygon = numel(domainMesh.connect);                                                     %Number of polygons
end

function domainMesh = read_mesh_wrench_domain_linelast2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the right circle boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.RightCircle = sscanf(line,'%d');
  % read indices of nodes that are located on the left half circle boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.LeftHalfCircle = sscanf(line,'%d');
  % read xmax, xmin, ymax, ymin for the bounding box
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs   
  % Right circle
  domainMesh.boundary_dofs.RightCircle=zeros(2*length(domainMesh.boundary_nodes.RightCircle),1);      
  range=1:length(domainMesh.boundary_nodes.RightCircle);
  domainMesh.boundary_dofs.RightCircle(2*range-1)=2*domainMesh.boundary_nodes.RightCircle-1;
  domainMesh.boundary_dofs.RightCircle(2*range)=2*domainMesh.boundary_nodes.RightCircle;    
  % Left half circle
  domainMesh.boundary_dofs.LeftHalfCircle=zeros(2*length(domainMesh.boundary_nodes.LeftHalfCircle),1);      
  range=1:length(domainMesh.boundary_nodes.LeftHalfCircle);
  domainMesh.boundary_dofs.LeftHalfCircle(2*range-1)=2*domainMesh.boundary_nodes.LeftHalfCircle-1;
  domainMesh.boundary_dofs.LeftHalfCircle(2*range)=2*domainMesh.boundary_nodes.LeftHalfCircle;       
end

function domainMesh = read_mesh_plate_with_hole_domain_linelast2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the lleft circle boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.Left = sscanf(line,'%d');
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.Right = sscanf(line,'%d');
  % read index of the node that is located on the middle of the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.MidNodeRightFace = sscanf(line,'%d');  
  % read xmax, xmin, ymax, ymin for the bounding box
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs   
  % Left
  domainMesh.boundary_dofs.Left=zeros(2*length(domainMesh.boundary_nodes.Left),1);      
  range=1:length(domainMesh.boundary_nodes.Left);
  domainMesh.boundary_dofs.Left(2*range-1)=2*domainMesh.boundary_nodes.Left-1;
  domainMesh.boundary_dofs.Left(2*range)=2*domainMesh.boundary_nodes.Left;    
  % Right
  domainMesh.boundary_dofs.Right=zeros(2*length(domainMesh.boundary_nodes.Right),1);      
  range=1:length(domainMesh.boundary_nodes.Right);
  domainMesh.boundary_dofs.Right(2*range-1)=2*domainMesh.boundary_nodes.Right-1;
  domainMesh.boundary_dofs.Right(2*range)=2*domainMesh.boundary_nodes.Right;
  % Node on the middle of the right boundary
  domainMesh.boundary_dofs.MidNodeRightFace=zeros(2*length(domainMesh.boundary_nodes.MidNodeRightFace),1);      
  range=1:length(domainMesh.boundary_nodes.MidNodeRightFace);
  domainMesh.boundary_dofs.MidNodeRightFace(2*range-1)=2*domainMesh.boundary_nodes.MidNodeRightFace-1;
  domainMesh.boundary_dofs.MidNodeRightFace(2*range)=2*domainMesh.boundary_nodes.MidNodeRightFace;  
end

function domainMesh = read_mesh_custom_linelast2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the dirichlet boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.Dirichlet = sscanf(line,'%d');
  % read indices of nodes that are located on the Neumann boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.Neumann = sscanf(line,'%d');
  % read xmax, xmin, ymax, ymin for the bounding box
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs   
  % Dirichlet boundary
  domainMesh.boundary_dofs.Dirichlet=zeros(2*length(domainMesh.boundary_nodes.Dirichlet),1);      
  range=1:length(domainMesh.boundary_nodes.Dirichlet);
  domainMesh.boundary_dofs.Dirichlet(2*range-1)=2*domainMesh.boundary_nodes.Dirichlet-1;
  domainMesh.boundary_dofs.Dirichlet(2*range)=2*domainMesh.boundary_nodes.Dirichlet;    
  % Neumann boundary
  domainMesh.boundary_dofs.Neumann=zeros(2*length(domainMesh.boundary_nodes.Neumann),1);      
  range=1:length(domainMesh.boundary_nodes.Neumann);
  domainMesh.boundary_dofs.Neumann(2*range-1)=2*domainMesh.boundary_nodes.Neumann-1;
  domainMesh.boundary_dofs.Neumann(2*range)=2*domainMesh.boundary_nodes.Neumann;       
end

function domainMesh = read_mesh_custom_without_boundary_data_linelast2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read xmax, xmin, ymax, ymin for the bounding box
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
end

function domainMesh = read_mesh_rectangular_domain_poisson2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');  
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the bottom boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.bottom = sscanf(line,'%d');
  % read indices of nodes that are located on the top boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.top = sscanf(line,'%d');
  % read indices of nodes that are located on the left boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.left = sscanf(line,'%d');    
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.right = sscanf(line,'%d');  
  % all boundary nodes
  domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.bottom;domainMesh.boundary_nodes.top;...
                             domainMesh.boundary_nodes.left;domainMesh.boundary_nodes.right];
  domainMesh.boundary_nodes.all=unique(domainMesh.boundary_nodes.all); % delete repeated nodes from boundary_nodes.all
  % read xmax, xmin, ymax, ymin for the rectangular domain
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs
  % all
  domainMesh.boundary_dofs.all=zeros(length(domainMesh.boundary_nodes.all),1);      
  range=1:length(domainMesh.boundary_nodes.all);
  domainMesh.boundary_dofs.all(range)=domainMesh.boundary_nodes.all;  
  % bottom
  domainMesh.boundary_dofs.bottom=zeros(length(domainMesh.boundary_nodes.bottom),1);      
  range=1:length(domainMesh.boundary_nodes.bottom);
  domainMesh.boundary_dofs.bottom(range)=domainMesh.boundary_nodes.bottom;   
  % top
  domainMesh.boundary_dofs.top=zeros(length(domainMesh.boundary_nodes.top),1);      
  range=1:length(domainMesh.boundary_nodes.top);
  domainMesh.boundary_dofs.top(range)=domainMesh.boundary_nodes.top;   
  % left
  domainMesh.boundary_dofs.left=zeros(length(domainMesh.boundary_nodes.left),1);      
  range=1:length(domainMesh.boundary_nodes.left);
  domainMesh.boundary_dofs.left(range)=domainMesh.boundary_nodes.left;  
  % right
  domainMesh.boundary_dofs.right=zeros(length(domainMesh.boundary_nodes.right),1);      
  range=1:length(domainMesh.boundary_nodes.right);
  domainMesh.boundary_dofs.right(range)=domainMesh.boundary_nodes.right;
  domainMesh.nvertex  = numel(domainMesh.coords) / 2;                                                  %Number of vertices
  domainMesh.npolygon = numel(domainMesh.connect);                                                     %Number of polygons
end

function domainMesh = read_mesh_wrench_domain_poisson2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the right circle boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.RightCircle = sscanf(line,'%d');
  % read indices of nodes that are located on the left half circle boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.LeftHalfCircle = sscanf(line,'%d');
  % read xmax, xmin, ymax, ymin for the bounding box
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs   
  % Right circle
  domainMesh.boundary_dofs.RightCircle=zeros(length(domainMesh.boundary_nodes.RightCircle),1);      
  range=1:length(domainMesh.boundary_nodes.RightCircle);
  domainMesh.boundary_dofs.RightCircle(range)=domainMesh.boundary_nodes.RightCircle;   
  % Left half circle
  domainMesh.boundary_dofs.LeftHalfCircle=zeros(length(domainMesh.boundary_nodes.LeftHalfCircle),1);      
  range=1:length(domainMesh.boundary_nodes.LeftHalfCircle);
  domainMesh.boundary_dofs.LeftHalfCircle(range)=domainMesh.boundary_nodes.LeftHalfCircle;
end

function domainMesh = read_mesh_plate_with_hole_domain_poisson2d(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the left boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.Left = sscanf(line,'%d');
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.Right = sscanf(line,'%d');
  % read index of the node that is located on the middle of the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.MidNodeRightFace = sscanf(line,'%d');  
  % read xmax, xmin, ymax, ymin for the bounding box
  line = fgets(mesh_file); % read commented line   
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs   
  % Left circle
  domainMesh.boundary_dofs.Left=zeros(length(domainMesh.boundary_nodes.Left),1);      
  range=1:length(domainMesh.boundary_nodes.Left);
  domainMesh.boundary_dofs.Left(range)=domainMesh.boundary_nodes.Left;   
  % Right
  domainMesh.boundary_dofs.Right=zeros(length(domainMesh.boundary_nodes.Right),1);      
  range=1:length(domainMesh.boundary_nodes.Right);
  domainMesh.boundary_dofs.Right(range)=domainMesh.boundary_nodes.Right;      
  % Node on the middle of the right boundary
  domainMesh.boundary_dofs.MidNodeRightFace=zeros(length(domainMesh.boundary_nodes.MidNodeRightFace),1);      
  range=1:length(domainMesh.boundary_nodes.MidNodeRightFace);
  domainMesh.boundary_dofs.MidNodeRightFace(range)=domainMesh.boundary_nodes.MidNodeRightFace;
end

end

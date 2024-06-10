%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_eig
%
% Created by : U. Zerinati 
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close, clc

fix_path();


%% PARAMETERS OF THE PDE
matProps.E      = 70;
matProps.nu     = 0.2;
matProps.mu     = matProps.E / (2 * (1 + matProps.nu));
matProps.lambda = matProps.E * matProps.nu / ( (1 + matProps.nu) * (1 - 2*matProps.nu));

fprintf('\n\nParameters: ')                                                                          %Print the parameters of the PDE 
fprintf('\nYung modulus =  %.2f\n',matProps.E)
fprintf('\nPoisson ratio =  %.2f\n',matProps.nu)
fprintf('\nLame mu =  %.2f\n',matProps.mu)
fprintf('\nLame lambda =  %.2f\n',matProps.lambda)

%% DEFINITION OF THE FUNCTIONS
[f1, f2, u1, u1_x, u1_y, u2, u2_x, u2_y] = problem_test_elasticity(4,matProps);                                   %Obtain problem functions

%% INFORMATION ON THE POLYNOMIALS
k = 1;                                                                                               %Degree of the polynomials used to solve the equation

fprintf('\nPolynomials degree for solving the equation: %d',k)

if (k ~= 1)
    error("Actually, the method only works with k = 1")
end

mesh = ["polygon_16.txt" "polygon_64.txt" "polygon_256.txt" "polygon_1024.txt" "polygon_4096.txt"];
%mesh = ["star_2_0.1.txt" "star_4_0.1.txt" "star_8_0.1.txt" "star_16_0.1.txt"];

%% INITIALIZE MEMORY
eigLightning = zeros(numel(mesh),10);
errLightning = zeros(numel(mesh),10);
rateLightning = zeros(numel(mesh),10);
diamVec = zeros(numel(mesh),1);

%% INFORMATION ON THE EIGENVALUE PROBLEM
eigExact = [1007.875293633368,...
            1007.875295841894,...
            1492.376174897711,...
            2165.996543952849,...
            2755.460012318369,...
            2882.723160229796,...
            2882.723160442838,...
            3529.696562815136,...
            4082.412001289926,...
            4082.412020260300];

numEigs = 10; meshIndex = 1;

for mesh_filename = mesh

tic

%% PRINT INIT MESSAGE TO SCREEN
fprintf('\n\n[%.2f] Starting the method... \n',toc);

%% READ THE MESH
fprintf('[%.2f] Reading mesh [%12s]...\n',toc, mesh_filename);
domainMesh    = read_mesh(mesh_filename);                                                            %Read mesh

%% OBTAIN INFORMATION ON THE MESH
fprintf('[%.2f] Obtaining information on the mesh... \n',toc);
                                                    
domainMesh          = add_edges(domainMesh, k);                                                   %Add to the domainMesh structure the edges
boundary_vertex     = domainMesh.boundary_nodes.all;                                                 %Extract boundary nodes

domainMesh.boundary_edges = get_boundary_edges(domainMesh);                                          %Get the boundary edges
domainMesh.internal_edges = setdiff(1:domainMesh.nedges, domainMesh.boundary_edges);                 %Get the internal edges

[boundary_dofs, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k);

%% ASSEMBLYING LIGHTNING VEM MATRICES
fprintf('[%.2f] Assemblying element matrices...\n',toc); 
[K_global, M_global, f_global, domainMesh] = elasticity_assembly(domainMesh, matProps, f1, f2, k);

boundary_dofs = [boundary_vertex; boundary_vertex + domainMesh.nvertex];
internal_dofs = setdiff(1:size(K_global,1), boundary_dofs);                                          %Get the indexes of the internal edges

diamVec(meshIndex) = domainMesh.diameter;

AII = K_global(internal_dofs,internal_dofs);                                                                %Matrix internal - internal
MII = M_global(internal_dofs,internal_dofs);

fprintf('[%.2f] Solving the eigenvalue problem...\n',toc); 

eigLightning(meshIndex,:) = eigs(AII,MII,numEigs,'smallestabs');
errLightning(meshIndex,:) = abs(eigLightning(meshIndex,:) - eigExact);

if meshIndex > 1
    for eigIndex = 1:numEigs
        rate = polyfit(log(diamVec([meshIndex-1,meshIndex])), log(errLightning([meshIndex-1,meshIndex],eigIndex)),1);
        rateLightning(meshIndex,eigIndex) = rate(1);
    end
end

meshIndex = meshIndex + 1;


end

nus = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.4999];
%% INITIALIZE MEMORY
eigPoisson = zeros(numel(mesh),numel(nus));
errPoisson = zeros(numel(mesh),numel(nus));
diamVec = zeros(numel(mesh),1);

%% INFORMATION ON THE EIGENVALUE PROBLEM
eigExact = [1007.875293633368, 1043.4500230981541, 1109.5419789873215, 1230.7475057513004,...
            1295.3480604990973, 1256.4697090150764, 1221.442537732282]

poissonIndex = 1; meshIndex = 1;

for nu = nus

    matProps.E      = 70;
    matProps.nu     = nu;
    matProps.mu     = matProps.E / (2 * (1 + matProps.nu));
    matProps.lambda = matProps.E * matProps.nu / ( (1 + matProps.nu) * (1 - 2*matProps.nu));
    
    for mesh_filename = mesh

        [f1, f2, u1, u1_x, u1_y, u2, u2_x, u2_y] = problem_test_elasticity(4,matProps);

        %% INFORMATION ON THE POLYNOMIALS
        k = 1;
        %% PRINT INIT MESSAGE TO SCREEN
        fprintf('\n\n[%.2f] Starting the method... \n',toc);

        %% READ THE MESH
        fprintf('[%.2f] Reading mesh [%12s]...\n',toc, mesh_filename);
        domainMesh    = read_mesh(mesh_filename);                                                            %Read mesh

        %% OBTAIN INFORMATION ON THE MESH
        fprintf('[%.2f] Obtaining information on the mesh... \n',toc);
                                                            
        domainMesh          = add_edges(domainMesh, k);                                                   %Add to the domainMesh structure the edges
        boundary_vertex     = domainMesh.boundary_nodes.all;                                                 %Extract boundary nodes

        domainMesh.boundary_edges = get_boundary_edges(domainMesh);                                          %Get the boundary edges
        domainMesh.internal_edges = setdiff(1:domainMesh.nedges, domainMesh.boundary_edges);                 %Get the internal edges

        [boundary_dofs, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k);

        %% ASSEMBLYING LIGHTNING VEM MATRICES
        fprintf('[%.2f] Assemblying element matrices...\n',toc); 
        [K_global, M_global, f_global, domainMesh] = elasticity_assembly(domainMesh, matProps, f1, f2, k);

        boundary_dofs = [boundary_vertex; boundary_vertex + domainMesh.nvertex];
        internal_dofs = setdiff(1:size(K_global,1), boundary_dofs);                                          %Get the indexes of the internal edges

        diamVec(meshIndex) = domainMesh.diameter;

        AII = K_global(internal_dofs,internal_dofs);                                                                %Matrix internal - internal
        MII = M_global(internal_dofs,internal_dofs);

        fprintf('[%.2f] Solving the eigenvalue problem...\n',toc); 

        eigPoisson(meshIndex,poissonIndex) = eigs(AII,MII,1,'smallestabs');
        errPoisson(meshIndex,poissonIndex) = abs(eigPoisson(meshIndex,poissonIndex) - eigExact(poissonIndex));
        meshIndex = meshIndex + 1;
    end
    poissonIndex = poissonIndex + 1;
end

loglog(diamVec,errPoisson(:,1), 'k*-','LineWidth',2, 'MarkerSize', 10)
hold on
loglog(diamVec,errPoisson(:,3), 'r*-','LineWidth',2, 'MarkerSize', 10)
loglog(diamVec,errPoisson(:,5),'b*-','LineWidth',2, 'MarkerSize', 10)
loglog(diamVec,errPoisson(:,7),'m*-','LineWidth',2, 'MarkerSize', 10)
grid on
axis square
legend('$\nu=0.2$', '$\nu=0.3$', '$\nu=0.4$', '$\nu=0.4999$', 'interpreter', 'latex')
xlabel('$h$','interpreter','latex')
ylabel('$|\lambda_1-\lambda_{1,h}|$','interpreter','latex')
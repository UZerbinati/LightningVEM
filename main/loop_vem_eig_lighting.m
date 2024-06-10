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
matProps.sigma   = 0;                                                                                %Reaction  coefficient
matProps.epsilon = 1;                                                                                %Diffusion coefficient
matProps.beta{1} = @(x,y) 0.*x + 0.*y;
matProps.beta{2} = @(x,y) 0.*x + 0.*y;

beta1      = func2str(matProps.beta{1});                                                             %Convert functions to string
beta2      = func2str(matProps.beta{2});
beta1(1:6) = [];                                                                                     %Remove @(x,y) 
beta2(1:6) = [];

fprintf('\n\nParameters: ')                                                                          %Print the parameters of the PDE 
fprintf('\nEpsilon =  %.2f\n',matProps.epsilon)
disp(['Beta    = [ ' beta1 ' ; ' beta2 ' ]'])
fprintf('Sigma   =  %.2f\n',matProps.sigma)

%% DEFINITION OF THE FUNCTIONS
[f, g, u, grad_u_x, grad_u_y] = problem_test_lighting(1,matProps);                                   %Obtain problem functions

%% INFORMATION ON THE POLYNOMIALS
k = 1;                                                                                               %Degree of the polynomials used to solve the equation

%polynomial = get_polynomial_info(k);
fprintf('\nPolynomials degree for solving the equation: %d',k)

if (k ~= 1)
    error("Actually, the method only works with k = 1")
end

%mesh = ["polygon_16.txt" "polygon_64.txt" "polygon_256.txt" "polygon_1024.txt"];
%mesh = ["star_2_0.1.txt" "star_4_0.1.txt" "star_8_0.1.txt" "star_16_0.1.txt"];
mesh = ["polygon_1024.txt"];

%% INFORMATION ON THE EIGENVALUE PROBLEM
eigExact = [2 5 5 8 10 10 13 13 17 17 18 20 20 25 25 26 26 29 29 32];
numEigs = 20; meshIndex = 1;

%% INITIALIZE MEMORY
eigLightning = zeros(numel(mesh),numEigs);
errLightning = zeros(numel(mesh),numEigs);
rateLightning = zeros(numel(mesh),numEigs);
diamVec = zeros(numel(mesh),1);


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
[K_global, M_global, f_global, domainMesh] = vem_lighting_assembly(domainMesh, matProps, f, k);

internal_dofs = setdiff(1:size(K_global,1), boundary_dofs);                                          %Get the indexes of the internal edges

diamVec(meshIndex) = domainMesh.diameter;

AII = K_global(internal_dofs,internal_dofs);                                                                %Matrix internal - internal
MII = M_global(internal_dofs,internal_dofs);

fprintf('[%.2f] Solving the eigenvalue problem...\n',toc); 

eigLightning(meshIndex,:) = eigs(AII,MII,numEigs,'smallestabs') ./ pi ./ pi;
errLightning(meshIndex,:) = abs(eigLightning(meshIndex,:) - eigExact);

if meshIndex > 1
    for eigIndex = 1:numEigs
        rate = polyfit(log(diamVec([meshIndex-1,meshIndex])), log(errLightning([meshIndex-1,meshIndex],eigIndex)),1);
        rateLightning(meshIndex,eigIndex) = rate(1);
    end
end

meshIndex = meshIndex + 1;


end

polynomial = get_polynomial_info(k);
%% ASSEMBLYING VEM MATRICES
fprintf('[%.2f] Assemblying element matrices...\n',toc); 
[A1, A2, B1, B2, f_global, domainMesh] = vem_adr_assembly(domainMesh, matProps, polynomial, f);

A1 = A1(internal_dofs,internal_dofs);
A2 = A2(internal_dofs,internal_dofs);
B1 = B1(internal_dofs,internal_dofs);
B2 = B2(internal_dofs,internal_dofs);

numTest = 51;
a = linspace(0.00,10,numTest);
beta = 5;

eigVEM = zeros(numEigs,numTest);
for i = 1:numTest

    
    A = A1 + a(i)*A2;
    B = B1 + beta*B2;

    eigVEM(:,i) = eigs(A,B,numEigs,'smallestabs') ./ pi ./ pi;
end

for i = 1:numEigs
    grid on
    plot(a,eigVEM(i,:),'o-', 'Linewidth', 1.5, 'MarkerSize', 5);
    hold on
    plot(a(end)+10/numTest,eigLightning(1,i),'k*')
    plot(a(end)+10/numTest,eigExact(i),'ro')
    ylim([0,40])
    xlabel('$\alpha$','interpreter','latex')
    ylabel('$\lambda_i$','interpreter','latex')
end
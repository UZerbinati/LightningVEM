% vem_lighting: This function computes the numerical solution of the PDE
% -eps * div(grad(u)) + beta * grad(u) + sigma * u = f
% using virtual element method (VEM).
%
% The basis functions of the virtual element space are computed using the lighting technique.
%
% The user can set the parameters of the PDE in lines 29-32.
%
% The test problem is set in line 43.
% 
% The mesh is selected in line 61. The meshes are constructed using VEMLAB 2.4  
% (https://camlab.cl/software/vemlab/) and the functions plot_mesh2d and
% read_mesh are from that software.

%% INITIALIZATION
clear; close; clc;

fix_path();

tic 

%% PRINT INITIAL MESSAGE
fprintf('[%.2f] Solution of the Advection - Diffusiom - Reaction problem with VEM ',toc);
fprintf('\n-eps * div(grad(u)) + beta * grad(u) + sigma * u = f')              

%% PARAMETERS OF THE PDE
matProps.sigma   = 0;                                                                                %Reaction  coefficient
matProps.epsilon = 1;                                                                                %Diffusion coefficient
matProps.beta{1} = @(x,y) 0 + 0.*x + 0.*y;
matProps.beta{2} = @(x,y) 0 + 0.*x + 0.*y;

beta1      = func2str(matProps.beta{1});                                                             %Convert functions to string
beta2      = func2str(matProps.beta{2});
beta1(1:6) = [];                                                                                     %Remove @(x,y) 
beta2(1:6) = [];

fprintf('\n\nParameters: ')                                                                          %Print the parameters of the PDE 
fprintf('\nEpsilon =  %.2f\n',matProps.epsilon)
disp(['Beta    = [ ' beta1 ' ; ' beta2 ' ]'])
fprintf('Sigma   =  %.2f\n',matProps.sigma)

%% DEFINITION OF THE FUNCTIONS
[f, g, u, grad_u_x, grad_u_y] = problem_test_lighting(4,matProps);                                   %Obtain problem functions

%% INFORMATION ON THE POLYNOMIALS
k = 3;                                                                                               %Degree of the polynomials used to solve the equation

%polynomial = get_polynomial_info(k);
fprintf('\nPolynomials degree for solving the equation: %d',k)

%% PRINT INIT MESSAGE TO SCREEN
fprintf('\n\n[%.2f] Starting the method... \n',toc);

%% READ THE MESH
fprintf('[%.2f] Reading a mesh...\n',toc);
mesh_filename = 'polygon_256.txt';  
domainMesh    = read_mesh(mesh_filename);                                                            %Read mesh

%% OBTAIN INFORMATION ON THE MESH
fprintf('[%.2f] Obtaining information on the mesh... \n',toc);
                                                    
domainMesh          = add_edges(domainMesh, k);                                                      %Add to the domainMesh structure the edges
boundary_vertex     = domainMesh.boundary_nodes.all;                                                 %Extract boundary nodes

domainMesh.boundary_edges = get_boundary_edges(domainMesh);                                          %Get the boundary edges
domainMesh.internal_edges = setdiff(1:domainMesh.nedges, domainMesh.boundary_edges);                 %Get the internal edges

[boundary_dofs, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k);

plot_mesh2d(domainMesh);                                                                             %Mesh plot
                                                                       
%% ASSEMBLYING DIFFUSION CONVECTION MATRIX
fprintf('[%.2f] Assemblying element matrices...\n',toc); 
[K_global, M_global, f_global, domainMesh] = vem_lighting_assembly(domainMesh, matProps, f, k);

%% CONVERTING SYSTEM MATRIX TO A SPARSE MATRIX
A = sparse(K_global);

%% CALCULATE BOUNDARY DATAS
fprintf('[%.2f] Enforcing Dirichlet boundary conditions...\n',toc);                                  %Impose boundary nodes

f_global(boundary_vertex)  = g(domainMesh.coords(boundary_vertex,1), ...                            
                               domainMesh.coords(boundary_vertex,2));
f_global(boundary_intdofs) = g(domainMesh.intcoords(boundary_intdofs-domainMesh.nvertex,1), ...
                                   domainMesh.intcoords(boundary_intdofs-domainMesh.nvertex,2));

%% SOLVE THE SYSTEM
fprintf('[%.2f] Solving system of linear equations...\n',toc);  

internal_dofs = setdiff(1:size(K_global,1), boundary_dofs);                                          %Get the indexes of the internal edges

AII = A(internal_dofs,internal_dofs);                                                                %Matrix internal - internal
AIB = A(internal_dofs, boundary_dofs);                                                               %Matrix internal - boundary
fI  = f_global(internal_dofs);                                                                       %Interior load term 
UB  = f_global(boundary_dofs);                                                                       %Force boundary values
UI  = AII \ (fI - AIB * UB);                                                                         %Solve system

MII = M_global(internal_dofs,internal_dofs);
U                = zeros(size(A,1),1);                                                               %Solution vector
U(internal_dofs) = UI;
U(boundary_dofs) = UB;

U_ex                           = zeros(domainMesh.nvertex + domainMesh.nedges * (k-1), 1);           %Exact solution
U_ex(1:domainMesh.nvertex)     = u(domainMesh.coords(:,1),domainMesh.coords(:,2));

%% ERRORS COMPUTATION
fprintf('[%.2f] Computing errors...\n',toc);
[errL2, errH1] = compute_errors_lighting(domainMesh, U, u, grad_u_x, grad_u_y, k);

fprintf("\n[%.2f] Errors computed: ",toc)
fprintf("\nL2 norm    : %f", errL2)
fprintf("\nH1 seminorm: %f\n", errH1)

% save("A.mat",'AII')
% save("B.mat",'MII')

%% PLOT
connect = domainMesh.connect;
coords  = domainMesh.coords;

figure(1)
for i = 1:domainMesh.npolygon
    hold on
    fill3( coords(connect{i},1), coords(connect{i},2), U(connect{i}), mean(U(connect{i})))
    axis square
end
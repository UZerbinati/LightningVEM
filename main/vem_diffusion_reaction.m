%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_diffusion_convection_cip
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function solves the diffusion-convection problem using virtual element methods and the 
% continuous interior penalty
% 
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May 29, 2022: implemented CIP for k > 1
% May  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION
clear; 
close;
clc;

addpath("..\")                                                                                       %Path to user defined functions
addpath("..\assembly\")
addpath("..\elements\")
addpath("..\errors\")
addpath("..\evalutation\")
addpath("..\matrices\")
addpath("..\mesh_files\")
addpath("..\preprocessing\")
addpath("..\quadrature\")
addpath("..\test\")
addpath("..\utility")

tic

%% PRINT INITIAL MESSAGE
fprintf('[%.2f] Solution of the Diffusion-Convection problem with Virtual Element Method ',toc);
fprintf('and Continuous Interior penalty');
fprintf('\n-eps * div(grad(u)) + beta * grad(u) = f')              

%% PARAMETERS OF THE PDE
matProps.sigma      = 1;                                                                                %Don't change it!
matProps.epsilon    = 1;

%Cip Parameter
matProps.gamma      = 0;                                                                             %CIP parameter
matProps.gamma_perp = 0.0;

fprintf('\n\nParameters: ')  
fprintf('\nEpsilon =  %.2f\n',matProps.epsilon)
fprintf("Gamma   =  %.3f\n",matProps.gamma);

%% OPTIONS FOR SOLVING THE PROBLEM

opt.stiffness     = "PiNabla";
opt.stabilization = "DofiDofi";

%% DEFINITION OF THE FUNCTIONS
[f, g, u, grad_u_x, grad_u_y] = problem_test(2,matProps);                                            %Obtain problem functions

%% INFORMATION ON THE POLYNOMIALS
k = 1;                                                                                               %Degree of the polynomials used to solve the equation

polynomial = get_polynomial_info(k);
fprintf('\nPolynomials degree for solving the equation: %d',k)
                   
%% PRINT INIT MESSAGE TO SCREEN
fprintf('\n\n[%.2f] Starting the method... \n',toc);

%% READ THE MESH
fprintf('[%.2f] Reading a mesh...\n',toc);

mesh_filename = 'trisquare64.txt';  
domainMesh    = read_mesh(mesh_filename);                                                            %Read mesh

%% OBTAIN INFORMATION ON THE MESH
fprintf('[%.2f] Obtaining information on the mesh... \n',toc);
                                                    
%domainMesh          = fix_boundary_points(domainMesh);
%domainMesh          = add_points(domainMesh,6);
domainMesh          = add_edges(domainMesh, k);                                                      %Add to the domainMesh structure the edges
boundary_vertex     = domainMesh.boundary_nodes.all;                                                 %Extract boundary nodes

domainMesh.boundary_edges = get_boundary_edges(domainMesh);                                          %Get the boundary edges
domainMesh.internal_edges = setdiff(1:domainMesh.nedges, domainMesh.boundary_edges);                 %Get the internal edges

[boundary_dofs, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k);

plot_mesh2d(domainMesh);                                                                             %Mesh plot

                                                                        
%% ASSEMBLYING DIFFUSION CONVECTION MATRIX
fprintf('[%.2f] Assemblying element matrices...\n',toc); 
[K_global, f_global, domainMesh] = vem_diffusion_reaction_assembly(domainMesh, matProps, ...
                                                                      polynomial, f, g, opt);

%% CONVERTING SYSTEM MATRIX TO A SPARSE MATRIX
A = K_global;

%% CALCULATE BOUNDARY DATAS

fprintf('[%.2f] Enforcing Dirichlet boundary conditions...\n',toc);                              %Impose boundary nodes
    
f_global(boundary_vertex)  = g(domainMesh.coords(boundary_vertex,1), ...                            
                                   domainMesh.coords(boundary_vertex,2));
f_global(boundary_intdofs) = g(domainMesh.intcoords(boundary_intdofs-domainMesh.nvertex,1), ...
                                   domainMesh.intcoords(boundary_intdofs-domainMesh.nvertex,2));

%% SOLVE THE SYSTEM
fprintf('[%.2f] Solving system of linear equations...\n',toc);  

internal_dofs = setdiff(1:size(K_global,1), boundary_dofs);                                      %Get the indexes of the internal edges
    
AII = A(internal_dofs,internal_dofs);                                                            %Matrix internal - internal
AIB = A(internal_dofs, boundary_dofs);                                                           %Matrix internal - boundary
fI  = f_global(internal_dofs);                                                                   %Interior load term 
UB  = f_global(boundary_dofs);                                                                   %Force boundary values
UI  = AII \ (fI - AIB * UB);                                                                     %Solve system
    
U                = zeros(size(A,1),1);                                                           %Solution vector
U(internal_dofs) = UI;
U(boundary_dofs) = UB;

%% EXACT SOLUTION ON THE NODES
U_ex                           = zeros(domainMesh.nvertex + domainMesh.nedges * (k-1), 1);           %Exact solution
U_ex(1:domainMesh.nvertex)     = u(domainMesh.coords(:,1),domainMesh.coords(:,2));
U_ex(domainMesh.nvertex+1:end) = u(domainMesh.intcoords(:,1),domainMesh.intcoords(:,2));


%% ERRORS COMPUTATION
fprintf('[%.2f] Computing errors...\n',toc);
[errL2, errH1, values] = compute_errors(domainMesh,polynomial,U,u, grad_u_x, grad_u_y);

errLinf = max(abs(U(1:domainMesh.nvertex + domainMesh.nedges * (k-1))-U_ex));

fprintf("\n[%.2f] Errors computed: ",toc)
fprintf("\nL2 norm    : %f", errL2)
fprintf("\nH1 seminorm: %f", errH1)
fprintf("\nLinf norm  : %f\n\n", errLinf)

toc

%% POST PROCESSING
connect = domainMesh.connect;
coords  = domainMesh.coords;

for i = 1:domainMesh.npolygon
    hold on
    fill3( coords(connect{i},1), coords(connect{i},2), U(connect{i}), values(i))
end

axis([0 1 0 1])
axis("square")
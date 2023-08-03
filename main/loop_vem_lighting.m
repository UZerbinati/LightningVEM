%% INITIALIZATION
clear; 
close;
clc;

addpath("../")                                                                                       %Path to user defined functions
addpath("../assembly/")
addpath("../chebfun/")
addpath("../elements/")
addpath("../errors/")
addpath("../evalutation/")
addpath("../matrices/")
addpath("../mesh_files/")
addpath("../preprocessing/")
addpath("../poly2D/")
addpath("../quadrature/")
addpath("../test/")
addpath("../utility")

tic

%% PRINT INITIAL MESSAGE
fprintf('[%.2f] Solution of the Diffusiom - Reaction problem with Virtual Element Method ',toc);
fprintf('\n-eps * div(grad(u)) + sigma * u = f')              

%% PARAMETERS OF THE PDE
matProps.sigma   = 1;                                                                                %Reaction  coefficient
matProps.epsilon = 1;                                                                                %Diffusion coefficient
matProps.tol     = 1e-4;                                                                             %Tolerance of Laplace
matProps.h       = 1e-7;                                                                             %Step for the Finite Difference
matProps.beta{1}    = @(x,y) -2*pi.*sin(pi.*(x + 2*y));
matProps.beta{2}    = @(x,y)    pi.*sin(pi.*(x + 2*y));


fprintf('\n\nParameters: ')  
fprintf('\nEpsilon =  %.2f', matProps.epsilon)
fprintf('\nSigma   =  %.2f', matProps.sigma)


%% DEFINITION OF THE FUNCTIONS
[f, g, u, grad_u_x, grad_u_y] = problem_test_lighting(3,matProps);                                   %Obtain problem functions

%% INFORMATION ON THE POLYNOMIALS
k = 1;                                                                                               %Degree of the polynomials used to solve the equation

polynomial = get_polynomial_info(k);
fprintf('\nPolynomials degree for solving the equation: %d',k)

%% INPUT OF THE MESH FILE NAME     
mesh = ["polygon_4.txt" "polygon_16.txt" "polygon_64.txt" "polygon_256.txt"]% "polygon_1024.txt"];
%mesh = ["square_4.txt" "square_16.txt" "square_64.txt" "square_256.txt" "square_1024.txt" "polygon_4096.txt" ];
%mesh = ["polygon_4.txt" "polygon_16.txt" "polygon_64.txt" "polygon_256.txt"];

%% INITIALIZE MEMORY
errL2 = zeros(numel(mesh),1);
errH1 = zeros(numel(mesh),1);
diam  = zeros(numel(mesh),1);

counter    = 1;
errL2vec   = zeros(numel(mesh),1);
errH1vec   = zeros(numel(mesh),1);
errLinfvec = zeros(numel(mesh),1);
diamvec    = zeros(numel(mesh),1);

for mesh_filename = mesh

%% PRINT INIT MESSAGE TO SCREEN
fprintf('\n[%.2f] Starting the method... \n',toc);

%% READ THE MESH
fprintf('[%.2f] Reading a mesh...\n',toc);

domainMesh    = read_mesh(mesh_filename);                                                            %Read mesh

%% OBTAIN INFORMATION ON THE MESH
fprintf('[%.2f] Obtaining information on the mesh... \n',toc);
                                                    
domainMesh          = add_edges(domainMesh, k);                                                      %Add to the domainMesh structure the edges
boundary_vertex     = domainMesh.boundary_nodes.all;                                                 %Extract boundary nodes

domainMesh.boundary_edges = get_boundary_edges(domainMesh);                                          %Get the boundary edges
domainMesh.internal_edges = setdiff(1:domainMesh.nedges, domainMesh.boundary_edges);                 %Get the internal edges

[boundary_dofs, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k);

%plot_mesh2d(domainMesh);                                                                             %Mesh plot
                                                                       
%% ASSEMBLYING DIFFUSION CONVECTION MATRIX
fprintf('[%.2f] Assemblying element matrices...\n',toc); 
[K_global, f_global, domainMesh] = vem_lighting_assembly(domainMesh, matProps, polynomial, f, k);

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

U                = zeros(size(A,1),1);                                                               %Solution vector
U(internal_dofs) = UI;
U(boundary_dofs) = UB;

U_ex                           = zeros(domainMesh.nvertex + domainMesh.nedges * (k-1), 1);           %Exact solution
U_ex(1:domainMesh.nvertex)     = u(domainMesh.coords(:,1),domainMesh.coords(:,2));
U_ex(domainMesh.nvertex+1:end) = u(domainMesh.intcoords(:,1),domainMesh.intcoords(:,2));

%% ERRORS COMPUTATION
fprintf('[%.2f] Computing errors...\n',toc);
[errL2, errH1] = compute_errors_lighting(domainMesh, U,u, grad_u_x, grad_u_y,k,matProps.h);

errLinf = max(abs(U(1:domainMesh.nvertex)-U_ex));

fprintf("[%.2f] Errors computed: ",toc)
fprintf("\nL2 norm    : %f", errL2)
fprintf("\nH1 seminorm: %f", errH1)
fprintf("\nLinf norm  : %f\n", errLinf)
%fprintf("\nCondition  : %f\n\n", condest(K_global))

toc

errL2vec(counter) = errL2;
errH1vec(counter) = errH1;
diamvec(counter) = domainMesh.diameter;

counter = counter + 1;

end

%% PLOT OF THE EOC

loglog(diamvec,errL2vec,'r*-','LineWidth',2, 'MarkerSize', 10)
hold on
loglog(diamvec,errH1vec,'b*-','LineWidth',2, 'MarkerSize', 10)
loglog(diamvec,diamvec,'k--','LineWidth',2, 'MarkerSize', 10)
loglog(diamvec,diamvec.^2/12,'k--','LineWidth',2, 'MarkerSize', 10)
txt = texlabel('$O(h)$');   
text(diamvec(end),diamvec(end)*0.7,txt,'interpreter','latex','HorizontalAlignment','center')
txt = texlabel('$O(h^2)$');   
text(diamvec(end),diamvec(end).^2/15*0.9,txt,'interpreter','latex','HorizontalAlignment','center')
grid on
axis square
legend("L2", "H1")
xlabel('$h$','interpreter','latex')
ylabel('error','interpreter','latex')
function [A1_loc, A2_loc, B1_loc, B2_loc, f_local, polygon] = vem_adr_element(domainMesh, matProps, polynomial, index, f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_adr_element
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function construct the local linear system associated to the discretization of the 
% advection-diffusion-reaction problem and the associated right-hand side.
%
% Input
% =====
% domainMesh : The struct that contains the information on the mesh that we are using
% matProps   : The struct that contains the parameters of the equation
% polynomial : The struct that contains the information on the polynomials that we are using
% index      : index of the polygon
% f          : The right-hand side of the linear system
%
% Output
% ======
% K_local : The local matrix of the linear system
% f_local : The local right-hand side
% polygon : The information on the polygon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXTRACT THE VERTICES
local_vertex = domainMesh.connect{index};                                                            %Indexes of the local vertices
local_vertex = [local_vertex; local_vertex(1)];                                                      %This is to simplify the implementation
vertex       = domainMesh.coords(local_vertex,:);                                                    %Coordinates of the vertices

%% GET THE INFORMATIONS OF THE POLYGON
polygon = get_polygon_info(vertex);                                                                  %Information on the polygon

polygon.vertex = vertex;                                                                             %Add the following fields
polygon.size   = polygon.nedges + polygon.nedges*(polynomial.k - 1) + polynomial.int;            

%% DEFINITION OF THE BASE FUNCTIONS AND THEIR GRADIENTS
base = @(x,y,k) (x - polygon.centroid(1)).^k(:,1) .* (y - polygon.centroid(2)).^k(:,2) ...           %Polynomial basis function
                ./ polygon.diameter.^(sum(k));                                                       %Gradient of the polynomial basis function
%grad = @(x,y,k) grad_fun(x,y,k,polygon);

%% GET QUADRATURE NODES AND WEIGHTS
p_quad = polygon_quadrature (polygon, polynomial.n);                                                 %Construct the quadrature on the polygon
b_quad = boundary_quadrature(polygon, polynomial.n);                                                 %Construct the quadrature on the boundary

%% POINTWISE VALUES OF THE BASIS FUNCTION
base_val_int   = base_evaluation_interior(p_quad, polynomial, polygon);
grad_val_int   = grad_evaluation_interior(base_val_int, p_quad, polynomial, polygon);

base_val_bound = base_evaluation_boundary(b_quad, polynomial, polygon);
grad_val_bound = grad_evaluation_boundary(base_val_bound, b_quad, polynomial, polygon);

%% CONSTRUCT THE PINABLA(K) PROJECTION

%CONSTRUCT THE MATRIX G AND GTILDE
G = vem_matrix_G(base_val_bound, grad_val_int, p_quad, polynomial, polygon);                         %Gradient of a polynomial times gradient of a polynomial

%CONSTRUCT THE MATRIX B
B = vem_matrix_B(grad_val_bound, polynomial, polygon);                                               %Gradient of a virtual times gradient of a polynomial

%PISTAR PROJECTION
Pistar  = G \ B;                                                                                     %Coefficients of the Pistar Projection

%% CONSTRUCT THE PI0(K) PROJECTION

% CONSTRUCT THE MATRIX H 
H = vem_matrix_H(base_val_int, p_quad, polynomial);                                                  %Polynomial times a polynomial
                
% EVALUATE THE PISTAR PROJECTION
Pistar_val_int = Pi_evaluate(Pistar, base_val_int);

% CONSTRUCT THE MATRIX C
C = vem_matrix_C(base_val_int, Pistar_val_int, p_quad, polynomial, polygon);

% CONSTRUCT THE MATRIX Pi0
Pi0 = H \ C;                                                                                         %Coefficient of the Pi0 projection

% Evaluate the Pi0 projection
Pi0_val_int = Pi_evaluate(Pi0, base_val_int);                                                        %evaluation of the projection 

%% CONSTRUCT STABILIZATION TERM
D = vem_matrix_D(base, polynomial, polygon, p_quad, base_val_int);                                   %DOFs matrix

Pi      = D * Pistar;
I       = eye(polygon.size);
Iminus  = I - Pi;                                                                                    %Stabilization term                   

%% CONSTRUCT THE STIFFNESS MATRIX

% ADD CONSISTENT-STIFFNESS
Gtilde      = G;                                                                                     %Gtilde is equal to G except for the first row 
Gtilde(1,:) = zeros(1, polynomial.dim);                                                              %that is equal to zero

A1_loc = Pistar' * Gtilde * Pistar;

A2_loc = Iminus' * Iminus;

%% CONSTRUCT THE MASS MATRIX
B1_loc = Pi0' * H * Pi0;

Pi      = D * Pi0;
I       = eye(polygon.size);
Iminus  = I - Pi;                                                                                    %Stabilization term                   

B2_loc = polygon.area*(Iminus' * Iminus);

%% CONSTRUCT f_local
f_local   = zeros(polygon.size,1);                                                                   %Memory allocation
f_val_int = f_evaluation_interior(f, p_quad);                                                        %evalutation of f on the quadrature nodes

for i=1:polygon.size
    
    integrand  = Pi0_val_int(:,:,:,i) .* f_val_int;
    f_local(i) = quadrature_2D(p_quad, integrand,"Evaluated");                                       %Construction of the local load term
    
end

%% STORE INFORMATION
polygon.G      = G;                                                                                  %Store information
polygon.B      = B;
polygon.H      = H;
polygon.C      = C;
polygon.Pistar = Pistar;
polygon.Pi0    = Pi0;
    
end  
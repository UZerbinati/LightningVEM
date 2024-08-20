function [K_local, M_loc, f_local, polygon] = elasticity_element_cr(vertex, matProps, f1, f2, k)

% vem_lighting_element: This function constructs the local matrix
% associated to the discretization of the Advection-Diffusion-Equation with
% the lighting virtual element method.
%
% Input parameters:
%   vertex: the structure that contains the information on the mesh
% matProps: the list of parameters for the PDE;
%        f: the right-hand side of the PDE.
%        k: the degree of the polynomials (actaully = 1)
%
% Output parameters:
% K_local: the global matrix;
% f_local: the right-hand side of the linear system;
% polygon: the information on the polygon.

%% EXTRACT THE VERTICES
vertex = [vertex; vertex(1,:)];

%% GET THE INFORMATIONS OF THE POLYGON
polygon = get_polygon_info(vertex);                                                                  %Information on the polygon

polygon.vertex = vertex;                                                                             %Add the following fields
polygon.size   = k * polygon.nedges + k*(k-1)/2;    

%% COMPUTE THE VIRTUAL BASIS FUNCTION
coords    = complex(vertex(:,1), vertex(:,2));
u_vem     = cell(polygon.size,1);                                                                    %Memory allocation

 lobatto   = lobatto_points_1d(k+1);

zerosBCs = cell(polygon.nedges,1);                                                                   %Zeros Boundary condition

for j = 1 : polygon.nedges 
    zerosBCs{j,1} = @(z) 0 + 0.*z;
end


for i = 4 : polygon.nedges

    arclength = @(z) abs(z - coords(i))/abs(coords(i+1)-coords(i)) * 2 - 1;                          %Map the edge into [-1;1]

    BCs = zerosBCs;

    polyf = polyfit(lobatto, [1, zeros(1,k)], k);
    BCs{i,1} = @(z) polyval(polyf,arclength(z));                                                     %ith-edge

    if i ~= 1                                                                                        %Previous edge   
        
        arclength2 = @(z) abs(z - coords(i-1))/abs(coords(i)-coords(i-1)) * 2 - 1;

        polyf      = polyfit(lobatto, [zeros(1,k), 1], k);
        BCs{i-1,1} = @(z) polyval(polyf,arclength2(z));

    else

        arclength2 = @(z) abs(z - coords(end-1))/abs(coords(end)-coords(end-1)) * 2 - 1;

        polyf      = polyfit(lobatto, [zeros(1,k), 1], k);
        BCs{end,1} = @(z) polyval(polyf,arclength2(z));        

    end                                                                                              

    u_vem{i} = lightningLaplace(coords(1:end-1),BCs,'tol',1e-6,'beta', ...
                                polygon.angles./pi,'noarnoldi','noplots');

end

u_vem{1} = @(z) 1 + 0.*z;
u_vem{2} = @(z) ( real(z) - polygon.centroid(1) ) / polygon.diameter;
u_vem{3} = @(z) ( imag(z) - polygon.centroid(2) ) / polygon.diameter;

%% GET QUADRATURE NODES AND WEIGHTS
p_quad = polygon_quadrature(polygon, k+1);                                                           %Construct the quadrature on the polygon

%% POINTWISE VALUES OF THE BASIS FUNCTION
val_int = vem_evaluation_interior(polygon.size, u_vem, p_quad);                                      %Evaluates the basis function in the quadrature nodes

%% CONSTRUCT THE MATRIX
Vertex = zeros(polygon.nedges, polygon.size);

for i = 1:polygon.nedges

    for j = 1:polygon.size

        Vertex(i,j) = u_vem{j}(coords(i,:)); 
    end

end

%% CONSTRUCT THE OTHER FUNCTIONS
b1 = [1; zeros(polygon.size-1, 1)];
b2 = [0; 1; zeros(polygon.size-2, 1)];
b3 = [0; 0; 1; zeros(polygon.size-3, 1)];

coeff = zeros(polygon.size,3); 
coeff(:,1) = Vertex \ b1;  
coeff(:,2) = Vertex \ b2;
coeff(:,3) = Vertex \ b3;

val_int = adjust_values(polygon.size, coeff, val_int);

%% APPROXIMATION OF THE DERIVATIVE WITH A FIRST ORDER FINITE DIFFERENCE
h = polygon.diameter * 1e-7;                                                                         %Size of the finite difference

%Evaluates the basis function near the quadrature points
aux        = p_quad;
aux.xi     = aux.xi + h;
val_int_xp = vem_evaluation_interior(polygon.size, u_vem, aux);
val_int_xp = adjust_values(polygon.size, coeff, val_int_xp);

aux.xi     = aux.xi - 2*h;
val_int_xm = vem_evaluation_interior(polygon.size, u_vem, aux);
val_int_xm = adjust_values(polygon.size, coeff, val_int_xm);

aux        = p_quad;
aux.eta    = aux.eta + h;
val_int_yp = vem_evaluation_interior(polygon.size, u_vem, aux);
val_int_yp = adjust_values(polygon.size, coeff, val_int_yp);

aux.eta    = aux.eta - 2*h;
val_int_ym = vem_evaluation_interior(polygon.size, u_vem, aux);
val_int_ym = adjust_values(polygon.size, coeff, val_int_ym);

%Approximate the derivative
grad_val_int_x = (val_int_xp - val_int_xm) / (2*h);
grad_val_int_y = (val_int_yp - val_int_ym) / (2*h);

%% APPROSIMATE THE DIFFUSIVE TERM
Gs_lighting = vem_matrix_lighting_Gs(grad_val_int_x, grad_val_int_y, p_quad);                          %Computing the stiffness matrix
L_lighting  = vem_matrix_lighting_L(grad_val_int_x, grad_val_int_y, p_quad);
                                                        
%% ASSEMBLY THE GLOBAL MATRIX
K_local = 2*matProps.mu * Gs_lighting + matProps.lambda * L_lighting;

%% ASSEMBLYH THE RIGHT-HAND SIDE
f_local   = zeros(2*polygon.size,1);                                                                   %Memory allocation
f1_val_int = f_evaluation_interior(f1, p_quad);
f2_val_int = f_evaluation_interior(f2, p_quad);%evaluation of f on the quadrature nodes

for i=1:polygon.size
    
    integrand  = val_int(:,:,:,i) .* f1_val_int;
    f_local(i) = quadrature_2D(p_quad, integrand,"Evaluated");                                       %Construction of the local load term
    
    integrand  = val_int(:,:,:,i) .* f2_val_int;
    f_local(i+polygon.size) = quadrature_2D(p_quad, integrand,"Evaluated");
end

%% STORE INFORMATION
polygon.basis = u_vem;                                                                               %Store the virtual functions

M_aux = vem_matrix_lighting_M(val_int, p_quad);

M_loc = [M_aux zeros(size(M_aux,1)); zeros(size(M_aux,1)) M_aux];
end  
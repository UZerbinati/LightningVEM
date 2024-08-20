function [K_local, M, f_local, polygon] = vem_lighting_element(vertex, matProps, f, k)

warning('off','MATLAB:nearlySingularMatrix')

% vem_lighting_element: This function constructs the local matrix
% associated to the discretization of the Advection-Diffusion-Equation with
% the lighting virtual element method.
%
% Input parameters:
%   vertex: the structure that contains the information on the mesh
% matProps: the list of parameters for the PDE;
%        f: the right-hand side of the PDE.
%        k: the degree of the polynomials (actaully = 1,2,3)
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

nedges  = polygon.nedges;

polygon.dofs                      = zeros(2,polygon.nedges*k);
polygon.dofs(:, 1:polygon.nedges) = vertex(1:end-1,:)';

Val = eye(k*nedges); 

if k > 1 

    for i = 1 : nedges

        param_x = @(t) (t+1) .* (vertex(i+1,1)-vertex(i,1)) ./ 2 + vertex(i,1);
        param_y = @(t) (t+1) .* (vertex(i+1,2)-vertex(i,2)) ./ 2 + vertex(i,2);

        lobatto = lobatto_points_1d(k+1);
        lobatto = lobatto(2:end-1);

        index = nedges + (k-1)*(i-1) + (1:k-1);
        polygon.dofs(1,index) = param_x(lobatto);
        polygon.dofs(2,index) = param_y(lobatto);

    end
end

%% POLYNOMIAL
n_poly  = (k + 1) * (k + 2) / 2;
deg     = zeros(n_poly,2);

for i=1:k                                                                                           %For each 1D degree
   
    deg((i*(i+1)/2 + 1:i*(i+1)/2 + i +1), 1) = (i:-1:0);                                            %Degrees of the polynomials with degree exactly equal to k
    deg((i*(i+1)/2 + 1:i*(i+1)/2 + i +1), 2) = (0:i);
    
end

poly = @(x,y,k) ((x - polygon.centroid(1)) / polygon.diameter).^(k(1)) ...
       .*       ((y - polygon.centroid(2)) / polygon.diameter).^(k(2));

%% COMPUTE THE VIRTUAL BASIS FUNCTION
coords    = complex(vertex(:,1), vertex(:,2));
u_vem     = cell(polygon.size,1);                                                                    %Memory allocation
ux_vem    = cell(polygon.size,1);
uy_vem    = cell(polygon.size,1);

lobatto   = lobatto_points_1d(k+1);

zerosBCs = cell(polygon.nedges,1);                                                                   %Zeros Boundary condition

%% GET QUADRATURE NODES AND WEIGHTS
p_quad = polygon_quadrature(polygon, k+10);                                                          %Construct the quadrature on the polygon

%% ZERO BOUNDARY CONDITIONS
for j = 1 : polygon.nedges 
    zerosBCs{j,1} = @(z) 0 + 0.*z;
end

for i = 1 : k*polygon.nedges                                                                         %Basis functions associated to the 1st and 2nd DoFs

    if i <= nedges

        i_edge = i;
        i_int  = 0;
 
    else

        i_edge = fix((i - nedges - 1) / (k - 1)) + 1;
        i_int  = mod((i - nedges - 1) , k - 1) + 1;

    end

    if (i_int == 0)
   
        arclength = @(z) abs(z - coords(i))/abs(coords(i+1)-coords(i)) * 2 - 1;                       %Map the edge into [-1;1]

        BCs = zerosBCs;
    
        polyf = polyfit(lobatto, [1, zeros(1,k)], k);
        BCs{i,1} = @(z) polyval(polyf,arclength(z));                                                  %ith-edge
    
        if i ~= 1                                                                                     %Previous edge   
            
            arclength2 = @(z) abs(z - coords(i-1))/abs(coords(i)-coords(i-1)) * 2 - 1;
    
            polyf      = polyfit(lobatto, [zeros(1,k), 1], k);
            BCs{i-1,1} = @(z) polyval(polyf,arclength2(z));
    
        else
    
            arclength2 = @(z) abs(z - coords(end-1))/abs(coords(end)-coords(end-1)) * 2 - 1;
    
            polyf      = polyfit(lobatto, [zeros(1,k), 1], k);
            BCs{end,1} = @(z) polyval(polyf,arclength2(z));        
    
        end                                                                                              
    
        [u_vem{i}, ux_vem{i}, uy_vem{i}] = laplace(coords(1:end-1),BCs,'tol',    1e-10,   'noplots');

    else 
    
        arclength = @(z) abs(z - coords(i_edge))/abs(coords(i_edge+1)-coords(i_edge)) * 2 - 1;
        BCs = zerosBCs;                                                                                 %Imposition of the BCs
           
        zeros_vec            = zeros(k+1, 1);
        zeros_vec(i_int + 1) = 1;

        polyf         = polyfit(lobatto, zeros_vec, k);
        BCs{i_edge,1} = @(z) polyval(polyf,arclength(z));

        [u_vem{i}, ux_vem{i}, uy_vem{i}] = laplace(coords(1:end-1),BCs,'tol',     1e-10,   'noplots');

    end

end

%% INTERIOR DOFS
if (k == 2)
   
    u_aux  = @(x,y) x.^2/2; 
    ux_aux = @(x,y) x;
    uy_aux = @(x,y) 0.*x + 0.*y;
    
    for i = 1:polygon.nedges
    
        BCs{i,1} = @(z) u_aux(real(z),imag(z));
    
    end
    
    [AAA, AAAx, AAAy] = laplace(coords(1:end-1),BCs,'tol',    1e-10,   'noplots');

    aux  = @(z) u_aux(real(z),imag(z)) - AAA(z);

    constant = quadrature_2D(p_quad, @(x,y) aux(x+1i.*y));

    u_vem{k*nedges + 1}  = @(z) aux(z) ./ (constant) .* polygon.area;
    ux_vem{k*nedges + 1} = @(z) (ux_aux(real(z),imag(z)) - AAAx(z)) ./ (constant) .* polygon.area;
    uy_vem{k*nedges + 1} = @(z) (uy_aux(real(z),imag(z)) - AAAy(z)) ./ (constant) .* polygon.area;

elseif (k == 3)

    u_aux  = @(x,y) (x - polygon.centroid(1)).^2 ./ polygon.diameter^2; 
    ux_aux = @(x,y) 2*(x - polygon.centroid(1)) ./ polygon.diameter^2;
    uy_aux = @(x,y) 0.*x + 0.*y;
    
    for i = 1:polygon.nedges
    
        BCs{i,1} = @(z) u_aux(real(z),imag(z));
    
    end
    
    [AAA, AAAx, AAAy] = laplace(coords(1:end-1),BCs,'tol',    1e-10,   'noplots');
    u_vem{k*nedges + 1}  = @(z) u_aux(real(z),imag(z)) - AAA(z);
    ux_vem{k*nedges + 1} = @(z) ux_aux(real(z),imag(z)) - AAAx(z);
    uy_vem{k*nedges + 1} = @(z) uy_aux(real(z),imag(z)) - AAAy(z);

    u_aux  = @(x,y) (x - polygon.centroid(1)).^3 ./ polygon.diameter.^3; 
    ux_aux = @(x,y) 3*(x - polygon.centroid(1)).^2 ./ polygon.diameter.^3;
    uy_aux = @(x,y) 0.*x + 0.*y;
    
    for i = 1:polygon.nedges
    
        BCs{i,1} = @(z) u_aux(real(z),imag(z));
    
    end
    
    [AAA, AAAx, AAAy] = laplace(coords(1:end-1),BCs,'tol',    1e-10,   'noplots');
    u_vem{k*nedges + 2}  = @(z) u_aux(real(z),imag(z)) - AAA(z);
    ux_vem{k*nedges + 2} = @(z) ux_aux(real(z),imag(z)) - AAAx(z);
    uy_vem{k*nedges + 2} = @(z) uy_aux(real(z),imag(z)) - AAAy(z);

    %% PLOT

    u_aux  = @(x,y) (y - polygon.centroid(2)).^3 ./ polygon.diameter.^3; 
    ux_aux = @(x,y) 0.*x + 0.*y;
    uy_aux = @(x,y) 3*(y - polygon.centroid(2)).^2 ./ polygon.diameter.^3;
    
    for i = 1:polygon.nedges
    
        BCs{i,1} = @(z) u_aux(real(z),imag(z));
    
    end
    
    [AAA, AAAx, AAAy] = laplace(coords(1:end-1),BCs,'tol',    1e-10,   'noplots');
    u_vem{k*nedges + 3}  = @(z) u_aux(real(z),imag(z)) - AAA(z);
    ux_vem{k*nedges + 3} = @(z) ux_aux(real(z),imag(z)) - AAAx(z);
    uy_vem{k*nedges + 3} = @(z) uy_aux(real(z),imag(z)) - AAAy(z);

    Val2 = zeros(3);

    for i = 1:3
        for j = 1:3

            integrand = @(z) u_vem{k*nedges + j}(z) .* poly(real(z),imag(z),deg(i,:));
            Val2(i,j)  = quadrature_2D(p_quad, @(x,y) integrand(x+1i.*y));
        end
    end

    Coeff2 = Val2 \ (polygon.area*eye(3));

end



%% ADJUST THE BASIS FUNCTIONS
dimension = k*nedges;


Coeff = eye(dimension);
%% POINTWISE VALUES OF THE BASIS FUNCTION
val_int  = vem_evaluation_interior(polygon.size, u_vem,  p_quad);                                      %Evaluates the basis function in the quadrature nodes
valx_int = vem_evaluation_interior(polygon.size, ux_vem, p_quad);
valy_int = vem_evaluation_interior(polygon.size, uy_vem, p_quad);

%Coeff2 = eye(3);
[val_int, valx_int, valy_int] = rescale_values(val_int, valx_int, valy_int, Coeff, dimension);

if k == 3
    range = k*nedges + (1:3);

    [val_int(:,:,:,range), valx_int(:,:,:,range), valy_int(:,:,:,range)] = rescale_values(val_int(:,:,:,range), valx_int(:,:,:,range), valy_int(:,:,:,range), Coeff2, 3);
end

%% APPROSIMATE THE DIFFUSIVE TERM
G_lighting = vem_matrix_lighting_G(valx_int, valy_int, p_quad);

%% APPROSIMATE THE REACTIVE TERM
M = vem_matrix_lighting_M(val_int, p_quad);                                                          %Computing the mass matrix

%% APPROSSIMATE THE ADVECTIVE TERM
b_val1 = f_evaluation_interior(matProps.beta{1}, p_quad);                                            %Beta aevaluation
b_val2 = f_evaluation_interior(matProps.beta{2}, p_quad);

T = vem_matrix_lighting_T(val_int, valx_int, valy_int, b_val1, b_val2, p_quad);

%% ASSEMBLY THE GLOBAL MATRIX
K_local = matProps.sigma * M + matProps.epsilon * G_lighting + T;

%% ASSEMBLYH THE RIGHT-HAND SIDE
f_local   = zeros(polygon.size,1);                                                                   %Memory allocation
f_val_int = f_evaluation_interior(f, p_quad);                                                        %evaluation of f on the quadrature nodes

for i=1:polygon.size
    
    integrand  = val_int(:,:,:,i) .* f_val_int;
    f_local(i) = quadrature_2D(p_quad, integrand,"Evaluated");                                       %Construction of the local load term
    
end

%% STORE INFORMATION
polygon.basis = u_vem;                                                                               %Store the virtual functions
polygon.basisX = ux_vem;
polygon.basisY = uy_vem;
polygon.coeff  = Coeff;
if k==3
polygon.coeff2  = Coeff2;
end
end
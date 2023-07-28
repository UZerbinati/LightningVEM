function [K_local, f_local, polygon] = vem_lighting_element(vertex, matProps, f, k)

%% EXTRACT THE VERTICES
vertex = [vertex; vertex(1,:)];

%% GET THE INFORMATIONS OF THE POLYGON
polygon = get_polygon_info(vertex);                                                                  %Information on the polygon

polygon.vertex = vertex;                                                                             %Add the following fields
polygon.size   = k * polygon.nedges + k*(k-1)/2;    
nedges = polygon.nedges;

%% COMPUTE THE VIRTUAL BASIS FUNCTION
coords    = complex(vertex(:,1), vertex(:,2));
u_vem     = cell(polygon.size,1);                                                                    %Memory allocation
f_vem     = cell(polygon.size,1);

lobatto   = lobatto_points_1d(k+1);

zerosBCs = cell(polygon.nedges,1);                                                                   %Zeros Boundary condition

for j = 1 : polygon.nedges 
    zerosBCs{j,1} = @(z) 0 + 0.*z;
end


for i = 1 : polygon.nedges

    % First type of DOFs: Vertices
    arclength = @(z) abs(z - coords(i))/abs(coords(i+1)-coords(i)) * 2 - 1;                          %Map the edge into [-1;1]

    BCs = zerosBCs;

    polyf = polyfit(lobatto, [1, zeros(1,k)], k);
    BCs{i,1} = @(z) polyval(polyf,arclength(z));                                                     %ith-edge

    if i ~= 1                                                                                        %Previous edge   
        
        arclength2 = @(z) abs(z - coords(i-1))/abs(coords(i)-coords(i-1)) * 2 - 1;

        polyf = polyfit(lobatto, [zeros(1,k), 1], k);
        BCs{i-1,1} = @(z) polyval(polyf,arclength2(z));

    else

        arclength2 = @(z) abs(z - coords(end-1))/abs(coords(end)-coords(end-1)) * 2 - 1;

        polyf = polyfit(lobatto, [zeros(1,k), 1], k);
        BCs{end,1} = @(z) polyval(polyf,arclength2(z));        

    end                                                                                              

    [u_vem{i}, ~, f_vem{i}, ~, ~, ~] = laplace(coords(1:end-1),BCs,'tol',matProps.tol,'noplots');    %Construct the basis function associated to the DOF

    %Second type of DOFs
    for j = 1:k-1                                                                                    
        
        index = nedges + (i-1)*(k-1) + j;                                                            %Local index of the DOF

        BCs = zerosBCs;                                                                              %Imposition of the BCs
           
        zeros_vec      = zeros(k+1, 1);
        zeros_vec(j+1) = 1;

        polyf    = polyfit(lobatto, zeros_vec, k);
        BCs{i,1} = @(z) polyval(polyf,arclength(z));

        [u_vem{index}, ~, f_vem{index}, ~, ~, ~] = laplace(coords(1:end-1),BCs,'tol',matProps.tol,'noplots');

    end


end

u_aux = @(x,y) x.^2/2 + y;

for i = 1:nedges

    BCs{i,1} = @(z) u_aux(real(z),imag(z));

end

% if (k>1)
%     AAA = laplace(coords(1:end-1),BCs,'tol',matProps.tol,'noplots');
%     [u_vem{k*nedges + 1}, ~, f_vem{k*nedges + 1}, ~, ~, ~] = @(z) u_aux(real(z),imag(z)) - AAA(z);
% end

%% GET QUADRATURE NODES AND WEIGHTS
p_quad = polygon_quadrature (polygon, k+1);                                                            %Construct the quadrature on the polygon

%% POINTWISE VALUES OF THE BASIS FUNCTION
vem_val_int = vem_evalutation_interior(polygon.size, u_vem, p_quad);

%% APPROXIMATION OF THE DERIVATIVE WITH A SECOND ORDER FINITE DIFFERENCE
h = polygon.diameter * matProps.h;

%+h
aux            = p_quad;
aux.xi         = aux.xi + h;
vem_val_int_xp = vem_evalutation_interior(polygon.size, u_vem, aux);
%vem_val_int_xp = vem_evalutation_interior(polygon.size, f_vem, aux);

aux.xi = aux.xi - 2*h;
vem_val_int_xm = vem_evalutation_interior(polygon.size, u_vem, aux);
%vem_val_int_xm = vem_evalutation_interior(polygon.size, f_vem, aux);

aux            = p_quad;
aux.eta        = aux.eta + h;
vem_val_int_yp = vem_evalutation_interior(polygon.size, u_vem, aux);
%vem_val_int_yp = vem_evalutation_interior(polygon.size, f_vem, aux);

aux.eta = aux.eta - 2*h;
vem_val_int_ym = vem_evalutation_interior(polygon.size, u_vem, aux);
%vem_val_int_ym = vem_evalutation_interior(polygon.size, f_vem, aux);

% %+2h
% aux            = p_quad;
% aux.xi         = aux.xi + 2*h;
% vem_val_int_2xp = vem_evalutation_interior(polygon.size, u_vem, aux);
% 
% aux.xi = aux.xi - 4*h;
% vem_val_int_2xm = vem_evalutation_interior(polygon.size, u_vem, aux);
% 
% aux            = p_quad;
% aux.eta        = aux.eta + 3*h;
% vem_val_int_2yp = vem_evalutation_interior(polygon.size, u_vem, aux);
% 
% aux.eta = aux.eta - 4*h;
% vem_val_int_2ym = vem_evalutation_interior(polygon.size, u_vem, aux);

%Approximate the derivative
grad_val_int_x = (vem_val_int_xp - vem_val_int_xm) / (2*h);
grad_val_int_y = (vem_val_int_yp - vem_val_int_ym) / (2*h);

%grad_val_int_x =  (imag(vem_val_int_yp + conj(vem_val_int))) / h;
%grad_val_int_y = -(imag(vem_val_int_xp + conj(vem_val_int))) / h;

% grad_val_int_x = (-vem_val_int_2xp + 8*vem_val_int_xp - 8*vem_val_int_xm + vem_val_int_2xm) / (12*h);
% grad_val_int_y = (-vem_val_int_2yp + 8*vem_val_int_yp - 8*vem_val_int_ym + vem_val_int_2ym) / (12*h);

G_lighting = vem_matrix_G_lighting(grad_val_int_x, grad_val_int_y, p_quad, polygon);                 %Computing the stiffness matrix

%% CONSTRUCT THE MASS MATRIX
M = vem_matrix_lighting_M(vem_val_int, p_quad);                                                      %Computing the mass matrix
K_local = matProps.sigma * M + matProps.epsilon * G_lighting;

%% CONSTRUCT THE CONVECTIVE MATRIX
b_val1_int = f_evalutation_interior(matProps.beta{1}, p_quad);                                       %Beta aevalutation
b_val2_int = f_evalutation_interior(matProps.beta{2}, p_quad);

T          = vem_matrix_lighting_T(vem_val_int, grad_val_int_x, grad_val_int_y, b_val1_int, b_val2_int, p_quad);
%% CONSTRUCT f_local
f_local   = zeros(polygon.size,1);                                                                   %Memory allocation
f_val_int = f_evalutation_interior(f, p_quad);                                                       %Evalutation of f on the quadrature nodes

for i=1:polygon.size
    
    integrand  = vem_val_int(:,:,:,i) .* f_val_int;
    f_local(i) = quadrature_2D(p_quad, integrand,"Evaluated");                                       %Construction of the local load term
    
end

K_local = K_local + T;

%% STORE INFORMATION
polygon.basis = u_vem;                                                                               %Store the virtual functions
%polygon.f     = f_vem;
end  
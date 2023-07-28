function [p_quad] = polygon_quadrature(polygon, n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: polygon_quadrature
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function constructs a quadrature that is exact of order 2n-1 for a polygon
%
% Input
% =====
% polygon : The polygon
% n       : Number of nodes for the 1D Gauss quadrature
%
% Output
% ======
% p_quad : A struct that contains the following fields
%        
%        bool  : An array where 1 indicates that the edges is in the quadrature
%        xi    : x-coordinates
%        eta   : y-coordinates
%        omega : The weights of the quadrature
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% OBTAIN THE INDEXES I
vertex = polygon.vertex;
nedges = polygon.nedges;
alpha  = polygon.centroid(1);
ni     = zeros(nedges,1);

I_omega_alpha = zeros(nedges, 1);                                                                    %Memory allocation for I_omega_alpha

for i = 1:nedges                                                                                     %Construction of I_omega_alpha
   
    if ( (vertex(i,2)~=vertex(i+1,2)) && ( (vertex(i,1)~=alpha) || (vertex(i+1,1)~=alpha) ) )

        I_omega_alpha(i) = 1;
        
        if (vertex(i,1) == vertex(i+1,1))
            ni(i) = n;
        else
            ni(i) = n + 1;
        end

    end
    
end

%% OBTAIN n AND ni
tau_n     = gauss_points_1d(n);                                                                      %n Gauss nodes
tau_n1    = gauss_points_1d(n+1);                                                                    %n+1 Gauss nodes
lambda_n  = gauss_weights_1d(n);                                                                     %n Gauss weights
lambda_n1 = gauss_weights_1d(n+1);                                                                   %n+1 Gauss weighs

%% CONSTRUCT ETA AND XI
p_quad.bool  = I_omega_alpha;
p_quad.xi    = zeros(nedges,n+1,n);                                                                  %Memory allocation
p_quad.eta   = zeros(nedges,n+1);
p_quad.omega = zeros(nedges,n+1,n);

for i =1:nedges
    
    if (I_omega_alpha(i) ~=0)
        
        x_i = @(t) ( vertex(i+1,1) - vertex(i,1) ) / 2 .*t ...
                 + ( vertex(i+1,1) + vertex(i,1) ) / 2;
        y_i = @(t) ( vertex(i+1,2) - vertex(i,2) ) / 2 .*t ...
                 + ( vertex(i+1,2) + vertex(i,2) ) / 2;
        
        if ( ni(i) == n)

            p_quad.eta(i,1:n) = y_i(tau_n);
            p_quad.eta(i,n+1) = p_quad.eta(i,1);
            xi_taun  = x_i(tau_n);
            xma = xi_taun - alpha;
            xpa = xi_taun + alpha;

            for k = 1:n
                
                p_quad.xi(i,1:n,k)  = xma * tau_n(k) + xpa;
                p_quad.omega(i,1:n,k) = (vertex(i+1,2)-vertex(i,2)) .* xma .*lambda_n.*lambda_n(k); 
                p_quad.xi(i,n+1,k) = p_quad.xi(i,1,k); 
                         
            end

        else

            p_quad.eta(i,:)   = y_i(tau_n1);
            xi_taun1 = x_i(tau_n1);
            xma = xi_taun1 - alpha;
            xpa = xi_taun1 + alpha;

    %        p_quad.xi(i,:,:) = (xi_taun1 - alpha);
            
            for k = 1:n

                p_quad.xi(i,:,k) = xma * tau_n(k) + xpa;

                p_quad.omega(i,:,k) = (vertex(i+1,2)-vertex(i,2)) .* xma .*lambda_n1.*lambda_n(k); 

            end

        end
        
    end

end

p_quad.xi = p_quad.xi / 2;
p_quad.omega = p_quad.omega / 4;
 function [polynomial] = get_polynomial_info(k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: get_polynomial_info
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function stores the information of the 2D polynomials that will be used in the code 
%
% Input
% =====
% k : The maximum degree of the polynomials
%
% Output
% ======
% polynomial : A struct that contains the information with the following fields:
%
%       dim     : The space dimension
%       int     : The number of interior DOFs of each polygon 
%       k       : The maximum degree of the 2D polynomials
%       n       : The number of quadrature nodes to compute exactly the integral of a polynomial 
%                 with degree k^2
%       points  : The Gauss quadrature nodes on (0,1)
%       val     : A matrix that contains the values of the basis functions associated to the DOFs 
%                 on the quadrature points
%       weights : The Gauss quadrature weights on (0,1)
%      
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    if ( (k ~= floor(k)) || (k<1) )                                                                  %Check if k is a positive integer
        error('ERROR! Invalid degree of polynomials');
    end
    
    polynomial.dim       = (k + 1) * (k + 2) / 2;                                                    %Number of 2D polynomials of degree k
    polynomial.dim_minus = (k^2 + k) / 2;
    polynomial.int       = (k^2 - k) / 2;                                                            %Number of interior DOFs of each polygon
    polynomial.k         = k;                                                                        %Degree of the polynomials
    polynomial.n         = k + 1;                                                                    %Number of nodes for the quadrature formula    
    polynomial.lap       = zeros(polynomial.dim, 2);

    polynomial.deg = zeros(polynomial.dim,2);                                                        
    polynomial.val = zeros(k,polynomial.n);

    polynomial.lobatto   = (lobatto_points_1d(polynomial.n) + 1) / 2;                                         
    
    values  = diag(ones(polynomial.k,1),-1);                                                         %Values of the polynomials (1 on a DOF, 0 on the others)
    values  = values(2:end,:);   
    points  = gauss_points_1d(polynomial.n);                                                         %Gauss quadrature nodes on the interval (-1,1)
    points  = (1 + points)/2;                                                                        %Shift the nodes on the interval (0,1)
    
    polynomial.weights = gauss_weights_1d(polynomial.n)/2;                                           %Weights of the Gauss quadrature on (0,1)
    polynomial.points  = points;                                                                     %Stores the quadrature nodes
    polynomial.val     = zeros(polynomial.n,polynomial.k);

    for i=1:k                                                                                        %For each 1D degree
       
        p = polyfit(polynomial.lobatto,values(i,:),polynomial.k);                                    %Polynomial that is 1 on the k-th node, 0 on the others 
        polynomial.val(:,i) = polyval(p,points)';                                                    %Evaluate p on the quadrature nodes

        polynomial.deg((i*(i+1)/2 + 1:i*(i+1)/2 + i +1), 1) = (i:-1:0);                              %Degrees of the polynomials with degree exactly equal to k
        polynomial.deg((i*(i+1)/2 + 1:i*(i+1)/2 + i +1), 2) = (0:i);
        
    end

    for d=4:polynomial.dim

        if polynomial.deg(d,1) >= 2
            polynomial.lap(d,1) = d - 2*sum(polynomial.deg(d,:)) + 1;
        end

        if polynomial.deg(d,2) >= 2
            polynomial.lap(d,2) = d - 2*sum(polynomial.deg(d,:)) - 1;
        end

    end
    
end
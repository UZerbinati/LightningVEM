function [derx, dery] = grad_fun(x,y,k,polygon)

if (k(1) ~= 0)

    derx =  k(:,1) .* (x - polygon.centroid(1)) .^ (k(:,1)-1) ...
                   .* (y - polygon.centroid(2)) .^  k(:,2) ...
                   ./ polygon.diameter .^ (sum(k));

else

    derx =  0 + 0.*x + 0.*y;

end

if (k(2) ~= 0)

    dery = k(:,2) .* (x - polygon.centroid(1)) .^  k(:,1) ...
                  .* (y - polygon.centroid(2)) .^ (k(:,2)-1) ...
                  ./ polygon.diameter .^ sum(k);

else
        
    dery = 0 + 0.*x + 0.*y;
    
end

end
function [Bp] = vem_matrix_Bp(beta_val1_bound, beta_val2_bound, base_val_bound, polynomial, polygon)

Bp  = zeros(polynomial.dim, polygon.size);

edges   = [polygon.edges(end); polygon.edges];
vnormal = [polygon.vnormal(end,:); polygon.vnormal];
k = polynomial.k;

weights = edges * polynomial.weights;

for j = 1:polygon.nedges

    index = polygon.nedges + 1+(j-1)*(k-1): polygon.nedges + j*(k-1);

    beta_inv = vnormal(j,:)    * [beta_val1_bound(j,end:-1:1); beta_val2_bound(j,end:-1:1)] ... 
             .* polynomial.val(:,1)';
    beta_str = vnormal(j+1,:)  * [beta_val1_bound(j+1,:); beta_val2_bound(j+1,:)] ... 
             .* polynomial.val(:,1)';
    beta_int = vnormal(j+1,:)  * [beta_val1_bound(j+1,:); beta_val2_bound(j+1,:)] ... 
             .* polynomial.val(:,2:end)';

    for i = 1:polynomial.dim 
        
        Bp(i,j)     =           beta_inv .* base_val_bound(j,end:-1:1,i) *  weights(j,:)';

        Bp(i,j)     = Bp(i,j) + beta_str .* base_val_bound(j+1,:,i)      * weights(j+1,:)';

        Bp(i,index) =           beta_int .* base_val_bound(j+1,:,i)      * weights(j+1,:)';
       
    end

end

end
function [val] = adjust_values(siz,coeff,val)

aux1 = coeff(1,1) * val(:,:,:,1);
aux2 = coeff(1,2) * val(:,:,:,1);
aux3 = coeff(1,3) * val(:,:,:,1);

for j = 2:siz
    aux1 = aux1 + coeff(j,1) * val(:,:,:,j); 
    aux2 = aux2 + coeff(j,2) * val(:,:,:,j);
    aux3 = aux3 + coeff(j,3) * val(:,:,:,j); 
end

val(:,:,:,1) = aux1;
val(:,:,:,2) = aux2;
val(:,:,:,3) = aux3;

end

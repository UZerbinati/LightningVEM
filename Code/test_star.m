clear;
clc;
close all;

n    = 8;
perc = 10;
meshname = 'meshes/star_' + string(n) + '_' + string(perc) + '.mat';

if(~isfile(meshname))
    vem_star_domain(n,perc);
end

load(meshname)
u = vem(meshname,@square_domain_rhs, ...
    @square_domain_boundary_condition);
s = size(elements,1);
figure
for i=1:s
    disp(i);
    I = elements{i};
    V = vertices(I, :);
    C = complex(V(:,1),V(:,2));
    nBC = size(C,1);
    for j=1:nBC
        if j+1 <= nBC
            BCs{j,1} = @(z) linearBC(z,C(j),C(j+1),u(I(j)),u(I(j+1)));
        else
            BCs{j,1} = @(z) linearBC(z,C(j),C(1),u(I(j)),u(I(1)));
        end
    end
    laplace(C,BCs,'tol',1e-6);
end
axis("square"), colorbar
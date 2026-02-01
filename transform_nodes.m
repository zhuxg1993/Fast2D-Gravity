function [n_node] = transform_nodes(node,mod)
nz = size(mod,1);
nx = size(mod,2);
n_node = zeros(nx*nz,4);

n = 1;
for j = 1 : nx
    for i = 1 : nz
        n_node(n,1) = node(i,j);
        n_node(n,2) = node(i,j+1);
        n_node(n,3) = node(i+1,j+1);
        n_node(n,4) = node(i+1,j);
        n = n + 1;
    end
end
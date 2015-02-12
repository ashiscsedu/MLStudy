function A=Affinity_SCC(coord, dims, K, sigma)

W=reshape(coord,2*size(coord,2),size(coord,3))';

A = scc(W,dims,K,sigma);

end
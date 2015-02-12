function A=Affinity_SLBF(coord, final_clusters,sigma)

W=reshape(coord,2*size(coord,2),size(coord,3))';

A=lbfsc(W,final_clusters,3,sigma); %SLBF-MS


end
function mse = computing_average_L2_error(data, dim, idx, ctr, dir)
% This ia originally a function in the SCC code by Guangliang Chen, see http://www.math.duke.edu/~glchen/scc.html
D = size(data,2);

K = max(idx);

if length(dim) == 1 && K > 1
    dim = dim*ones(K,1);
end

if nargin<4
    [ctr,dir,dim]= computing_centers_and_bases(data,idx,dim);
end

mse = 0;
for k = 1:K
    cls_k = data((idx==k),:);
    n_k = size(cls_k,1);
    if n_k > dim(k)
        mse = mse + sum(sum(((cls_k - repmat(ctr{k,1},n_k,1))*(eye(D) - dir{k,1}'*dir{k,1})).^2,2));
    %else
    %    mse = Inf;
    end
end

mse = sqrt(mse/sum(idx>0));
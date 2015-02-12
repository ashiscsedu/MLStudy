function [ctr,dir,dims]= computing_centers_and_bases(data,inds,dims)
% This ia originally a function in the SCC code by Guangliang Chen, see http://www.math.duke.edu/~glchen/scc.html
%D = size(data,2);

K = max(inds);

if length(dims) == 1 && K>1
    dims = dims*ones(K,1);
end

% intialization
ctr = cell(K,1);
dir = cell(K,1);

for k = 1:K
    
    cls_k = data((inds==k),:);
    n_k = size(cls_k,1);

    if n_k >= dims(k)+1
       
       ctr{k,1} = mean(cls_k,1);  
       cls_k = cls_k - repmat(ctr{k,1},n_k,1);
       [uk,sk,vk] = svd(cls_k,0);
    
       %dims(k) = extract_true_dimension(diag(sk(1:dims(k),1:dims(k))).^2,.9);
    
       dir{k,1} = vk(:,1:dims(k))';
    
    end

end
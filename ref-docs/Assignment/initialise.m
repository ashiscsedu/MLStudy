
function [clusters,clust_size]=initialise(coord,t,clust_size,Npoints,Ntotal,type,initI)
% Initialise only from the t frame of the sequence

clusters={};
if exist('initI','var')
   %%pre-initialised by other method
   O={};
   for i=1:max(initI)
       O{i}=find(initI==i);
       
       v=size(O{i},1)/Ntotal; %the number of points in current cluster as percentage of the total
       K=floor(clust_size*v); %the number of Npoints-clusters  from currrent cluster as percentage of the total
       K=min(K,size(O{i},1));

        sindx=randsample(O{i},K); %random sample K centers from current cluster
        C=squeeze(coord(1:2,t,sindx))'; %get their coordinates
        D=pdist2(C,squeeze(coord(1:2,t,O{i}))')'; %get the distances between every point in the cluster and the centers
       for c=1:K
        [Dist, DI] =sort(D(:,c));
        clusters{end+1}=O{i}(DI(1:Npoints));
        end

   end
   clust_size=size(clusters,2);
   
   

   
   
else
    
    if strcmp(type,'kmeans')
        
        %% kmeans initialisation
        kmeans_type='sample'; %choices 'cluster', 'uniform', 'sample'
        min_KMclust_size=1;     %minimum number of points in the K-means clusters
        
        
        %kmeans clustering for initialisation
        [memb,C,~, D] = kmeans(squeeze(coord(1:2,t,:))',clust_size,'EmptyAction','drop','Replicates',10, 'Start', kmeans_type);
        
        
        %discard clusters with very few points
        idx=[];
        for i=1:clust_size
            if sum(memb==i) < min_KMclust_size
                idx(end+1)=i;
            end
        end
        D(:,idx)=[];
        C(idx,:)=[];
        clust_size=clust_size-length(idx);
        
        
        
        
    else
        %% random sample initialisation
        sindx=randsample(Ntotal,clust_size);
        C=squeeze(coord(1:2,t,sindx))';
        D=pdist2(C,squeeze(coord(1:2,t,:))')';
        
        
    end
    
    
    %take Npoints-closer points to each cluster center
    clusters=cell(1,clust_size); %clusters
    for i=1:clust_size
        [Dist, DI] =sort(D(:,i));
        clusters{i}=DI(1:Npoints);
    end
end

function A=Affinity_LCV(Npoints,clust_size,coord,sigma)


%sigma is the parameter of the kernel


t=1; %the frame where to do initialisation
Ntotal=size(coord,3);
frames=size(coord,2);
beta=1;

%setup the basis views. Fixed for all data
basisView1_frame=1; %first frame
basisView2_frame=frames;%last frame




    if clust_size>Ntotal, clust_size=Ntotal; end
    if clust_size<=40
        [clusters,clust_size]=initialise(coord,t,clust_size,Npoints,Ntotal,'kmeans');
    else
        [clusters,clust_size]=initialise(coord,t,clust_size,Npoints,Ntotal,'other');
    end


%% Calculate LCV coefficients
Hubertau=15;


%generate a "motion hypothesis" for each cluster. In other words calculate
%the LCV coefficients for each cluster
a0=zeros(frames,clust_size); a1=zeros(frames,clust_size); a2=zeros(frames,clust_size); a3=zeros(frames,clust_size); a4=zeros(frames,clust_size);
b0=zeros(frames,clust_size); b1=zeros(frames,clust_size); b2=zeros(frames,clust_size); b3=zeros(frames,clust_size); b4=zeros(frames,clust_size);

for i=1:clust_size
    

    clusterBV1=[reshape(coord(1,basisView1_frame,clusters{i}),Npoints,1)  reshape(coord(2,basisView1_frame,clusters{i}),Npoints,1)];
    clusterBV2=[reshape(coord(1,basisView2_frame,clusters{i}),Npoints,1)  reshape(coord(2,basisView2_frame,clusters{i}),Npoints,1)];
    
    for f=1:frames %loop for all frames
        
        novelView_frame=f;
        
        
        clusterNV=[reshape(coord(1,novelView_frame,clusters{i}),Npoints,1) reshape(coord(2,novelView_frame,clusters{i}),Npoints,1)];
        
        S=[ones(size(clusterNV,1),1) clusterBV1  clusterBV2];
        C=S\clusterNV; %linear least squares
        
        a0(f,i)=C(1,1); a1(f,i)=C(2,1); a2(f,i)=C(3,1); a3(f,i)=C(4,1); a4(f,i)=C(5,1);
        b0(f,i)=C(1,2); b1(f,i)=C(2,2); b2(f,i)=C(3,2); b3(f,i)=C(4,2); b4(f,i)=C(5,2);
    end
    
    
end


%% Now classify each remaining point
%Use the pre-calculated coefficients to determine the best trajectory of
%each point in the scene

err=zeros(Ntotal,clust_size);
for j=1:Ntotal
    
    %the real point trajectories
    gt=[coord(1,:,j)' coord(2,:,j)'];
    
    %basis views
    objectBV1=gt(basisView1_frame,:);
    objectBV2=gt(basisView2_frame,:);
    
    
    %now loop for all the clusters and try to reconstruct the point
    %trajectories
    for i=1:size(clusters,2)
        %loop for all the frames
        novel=zeros(frames,2);
        for f=1:frames;
            novel(f,1)=a0(f,i)+a1(f,i)*objectBV1(1)+a2(f,i)*objectBV1(2)+a3(f,i)*objectBV2(1)+a4(f,i)*objectBV2(2);
            novel(f,2)=b0(f,i)+b1(f,i)*objectBV1(1)+b2(f,i)*objectBV1(2)+b3(f,i)*objectBV2(1)+b4(f,i)*objectBV2(2);
        end
        
        %%calculate the error for each cluster-reconstructed trajectory
        err(j,i)= (sqrt(1+norm(novel-gt,'fro')^2/Hubertau^2)-1)/frames;
        
    end
end





    KA= 1./sqrt(err.^2+sigma^2); %inverse multiquadric kernel.
    
    A=(KA*KA').^beta;
    A(1:Ntotal+1:Ntotal*Ntotal)=0; %enforce diagonal to 0
    
  



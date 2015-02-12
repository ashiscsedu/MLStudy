function [centers mincent ene]=rshlmnfnnbi(B,nc,cc,ss,opt,maxp,maxiter);

%[centers mincent ene]=rshlmnfnnbi(B,nc,cc,ss,maxp,maxiter);
% B: the data matrix of size N*D, representing N points in D-dimensional space

%nc number of centers

%cc is dimension

% ss is 0 for affine subspace or 1 for linear subspace

% opt==1 for LBF, and 2 for LBF-MS

%maxp is number of guess planes

%maxiter is maximum tries at finding best l1 subset of guess planes

%ene: the l1 energy

%centers: guessed planes

%mincent: output labels



nccc=ones(1,nc)*cc;% a vector representing the dimension of subspaces

if nargin<7
 %default value for maxp and maxiter   
    
    maxp=70*nc;
    
    maxiter=5*nc;
    
end

[N dim]=size(B);

maxp=min(N,maxp);% maxp can not be bigger than the number of points


rid=repmat(randpermute(floor(N)),1,1);

nB=sum(B.^2,2);

step=2;% find the best neighborhood by increasing its size by 2 each time

scales=40;%check 40 different neighborhoods at most
i=1;

for k=1:maxp% find 'maxp' candidates on subspaces
    % fit a good subspace according to option
        if opt==2
            [uid q]=gnhbd(B,rid(k),cc,2*cc,step,scales,ss,1,nB);
        elseif opt==1
            [uid q]=gnhbd2(B,rid(k),cc,2*cc,step,scales,ss,1,nB);
        elseif opt==3
            [uid q]=msgnhbd2(B,rid(k),cc,2*cc,step,scales,ss,1,nB);
        elseif opt==4
            [uid q]=msgnhbd(B,rid(k),cc,2*cc,step,scales,ss,1,nB);
        end
        
        
        scenters{i}=q;% the model of the i-th subspace. This is a matrix, the first row is the center, other rows are the basis of the subspace
        scentersdim(i)=cc;% the dimension of the subspace
        
        
        
        i=i+1;
    
    
    
end




    ncr=length(scenters);
    
    
    for k=1:ncr
        
        if ss==1
            
            ds(:,k)=kmccpdist(scenters{k},B,1);
            
        elseif ss==0
            
            ds(:,k)=kmccpdist(scenters{k},B);
        elseif ss==2
            ds(:,k)=gdistance(B,scenters{k})';
            
        end
        
    end
    inlierbound=0.1;
    
    ds=sqrt(ds);
    
    
    
    
    
    ene=inf;% the initialization of the l1 energy
    clear kk
    kk=find(scentersdim==cc(end));
    
    s=kk(ceil(rand*length(kk)));% pick any subspace candidate 
    
    used=s;
    
    for j=1:nc-1
        
        kk=find(scentersdim==nccc(end-j));

        nused=setdiff(kk,used);
        
        dm=min(ds(:,used),[],2);
        
        clear td
        
        for k=1:length(nused);
            
            td(k)=sum(min([dm ds(:,nused(k))],[],2));% calculate the l1 energy if the k-th unused subspace is added to the model
            
            
            
        end
        
        
        [ju juu]=min(td);
        used=cat(1,used,nused(juu));% pick the (j+1)-th subspace to minimize the l1 energy
        
    end
    
    nused=setdiff([1:ncr],used);%unused subspaces

    for it=1:maxiter
        
        tused=used;
        
        killid=ceil(rand*nc);
        
        tused(killid)=[];% remove one random subspace from the current model
        
        dm=min(ds(:,tused),[],2);
        
        clear td
        
        nused=setdiff((find(scentersdim==nccc(killid))),used);
        
        for k=1:length(nused);
            
            td(k)=sum(min([dm ds(:,nused(k))],[],2));% calculate the l1 energy if the k-th unused subspace is added to the model
            
            
            
        end
        
        
        
        [ju juu]=min(td);
        
        %ju
        
        if ju<ene% replace it by another subspace if it reduce the l1 energy
            
            used(killid)=nused(juu);
            
            ene=ju;
            
        end
        
        nused=setdiff([1:ncr],used);
        
        
    end
    
    

for k=1:nc
    
    centers{k}=scenters{used(k)};
    
end

for k=1:nc
    
    if ss==1
        
        d(:,k)=kmccpdist(centers{k},B,1);
        
    elseif ss==0
        
        d(:,k)=kmccpdist(centers{k},B);
        
        
    end
    
end
d=d.^0.5;%distances from points to the subspaces



[tmu mincent]=min(d,[],2);
ene=sum(tmu);% the final l1 energy

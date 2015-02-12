function fullDatabaseRun


rng('default');
addpath scc;
addpath rshlmnfnnbi;



path2seq='./Local_Hopkins155_slim';
DStruct=dir(path2seq);
Errors=[]; Times=[];
counter=1;
%loop for all the directories that contain the file "parameters.m"
for i=3:size(DStruct,1) %ignore . and ..
    
    
    if DStruct(i).isdir
        
        seqpath = fullfile(path2seq,DStruct(i).name);
        
        paramcheck=dir(fullfile(seqpath,'/parameters.m'));
        if size(paramcheck,1)~=0 %valid directory found
            
            
            %load data
            clear m snew complt;
            load([path2seq '/' DStruct(i).name '/Trajectories.mat']);
            load([path2seq '/' DStruct(i).name '/',DStruct(i).name,'_truth.mat']);
            run([path2seq '/' DStruct(i).name '/parameters.m']);
            
            
            
            %Do some processing about missing trajectories and outliers here
            
            %ignore outliers
            snew=s; snew(snew==0)=[];
            coord=coord(:,:,s>0);
            %ignore missing trajectories
            if exist('m','var')
                m(s==0,:)=[];
                complt=sum(m,2);
                snew(complt<frames)=[];
                coord=coord(:,:,complt==frames);
            else
                m=[];
            end
            s=snew;
            
            
            
            %%run the sequence
            fprintf('Running sequence %i: %s .... \n',counter, DStruct(i).name);
            
            
            
            
            %% EDIT FROM HERE ONWARDS
            
            
            
            
            % Affinity matrix 1
            p1=4; %the kernel parameter
            A1=Affinity_TK(coord,p1);
            
            % Affinity matrix 2
            p2=10; %the kernel parameter
            A2=Affinity_LCV(7,100*final_clusters,coord,p2);
            
            
            
            % Affinity matrix 3
            p3=0.01; %the kernel parameter
            
            A3=Affinity_SCC(coord,4,final_clusters,p3);
         
            
            
            % Affinity matrix 4
            p4=0.1; %the kernel parameter
            A4=Affinity_SLBF(coord,final_clusters,p4);
            
            
            
            
            % YOUR SPECTRAL CLUSTERING CODE HERE
            
            
         
            
            
            
            
            
            
            %Evaluation
            I=randi(final_clusters, size(coord,3),1); %random label vector. Replace this with your own clustering labels
            classification_error(s,I) %evaluate against the ground truth segmentation. Results in percentage

   
            counter=counter+1;
            
        end
        
    end
end




end







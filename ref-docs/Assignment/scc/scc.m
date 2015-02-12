function Aout = scc(X,d,K,sigma)



%global ABSOLUTE_MINIMUM
ABSOLUTE_MINIMUM = 1e-15;

OPTIONS = struct();
OPTIONS.n = d+2;
OPTIONS.c = K*100;

    if mod(OPTIONS.c,K)>0
        OPTIONS.c = K*ceil(OPTIONS.c/K);
        warning('The value of the parameter c in OPTIONS has been modified to be an integer multiple of K!'); %#ok<WNTAG>
    end


    OPTIONS.findOptimalSigma = 0;
    OPTIONS.sigma=sigma;

    OPTIONS.normalizeW = 1;


    OPTIONS.normalizeU = 1;


    OPTIONS.seedType = 'hard';



    OPTIONS.alpha = 0;








%% initialization

N = size(X,1);



% auxiliary 
averageL2Error = Inf;
sampleLabels = zeros(N,1);

if ~isfield(OPTIONS, 'initialLabels')
    % consider given data as one cluster
    sampleLabels1 = ones(N,1);
else
    % only for improving clusters obtained by other algorithms
    sampleLabels1 = OPTIONS.initialLabels;
end

averageL2Error1 = computing_average_L2_error(X, d, sampleLabels1);

%% main body of the code

q = OPTIONS.n; % q = d+2 by default
while averageL2Error1 < averageL2Error * 0.99 || q>2
%for iii = 1    
    sampleLabels = sampleLabels1;
    averageL2Error = averageL2Error1;
    
    %if isfield(OPTIONS,'sampledColumns')
    %    sampledColumns = OPTIONS.sampledColumns;
    %    OPTIONS = rmfield(OPTIONS,'sampledColumns');
    %else
    %   sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);
    %end

    sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);
    
    polarCurv = computing_polar_curvatures(X,sampledColumns,d);
    polarCurv_sorted = sort(reshape(polarCurv,1,N*OPTIONS.c));
    
    q = max(q-1,1);
    
        
        if isfield(OPTIONS,'sigma')
            sigma = OPTIONS.sigma;
        else
            sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K)),ABSOLUTE_MINIMUM);
            %sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K^(OPTIONS.n-1))),ABSOLUTE_MINIMUM);
        end

        isigma = 1/(2*sigma);
        A = exp(-polarCurv*isigma);
        sampleLabels1 = processing_affinities(A,K,OPTIONS);
        averageL2Error1 = computing_average_L2_error(X, d*ones(K,1), sampleLabels1);
        Aout=(A*A').^2;
        Aout(1:N+1:N*N)=0;
  
        
end % while averageL2Error1 < averageL2Error * 0.99

if averageL2Error1 < averageL2Error 
    averageL2Error = averageL2Error1;
    sampleLabels = sampleLabels1;
end

%do_plot_data(X,sampleLabels); title('clusters obtained by SCC');
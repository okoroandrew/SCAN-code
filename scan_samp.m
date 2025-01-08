
function [Omega,pred] = scan_samp(X,H,HatomInfo,Phi,threshold,rho,lambda,beta)

    [t,f] = size(X);
    
    % Draw random samples samples for stopping criteria
    size_ = 3;
    size_dist = 200000;
    r_end = t-size_+1;
    c_end = f-size_+1;
    distribution = zeros(1, size_dist);
    
    % Generate random indices for the top-left corner of the 5x5 matrices
    rowIndices = randi(r_end, [1, size_dist]);
    colIndices = randi(c_end, [1, size_dist]);
    
    for i = 1:size_dist
        subMatrice = X(rowIndices(i):(rowIndices(i) + size_-1), ...
                          colIndices(i):(colIndices(i) + size_-1));
        distribution(i) = mean(subMatrice(:));
    end
    
    pd = fitdist(distribution', 'Normal');
    
    % intilization
    Omega = ones(size(X));
    
    [p,f]= size(H);
    
    sum_H = sum(H,2);
    flat_X = sum(X);
    
    % subselecting H atom to speed up SCAN
    a =  (flat_X*H')./sum_H';
    [~,potFreqInd] = maxk(a,round(p*beta));
    H = H(potFreqInd, :);
    HatomInfo= HatomInfo(potFreqInd, :);
    
    % normalized H 
    Norm_H = vecnorm(H');
    
    pred{1} = zeros(t, f);
    
    i = 1;
    while true
    
        mu_X = sum(X.*Omega,2)./sum(Omega,2); 
        
        Z = (X-mu_X).*Omega;
    
        A =  Z*H'  ;
    
        cc = A./Norm_H;
        a = mean(cc);

        [~,m] = max(a);
    
        chosen_atom = H(m,:);
        atom_start = HatomInfo(m,2);
        atom_end = HatomInfo(m,3);
        
        ident_f = atom_start:atom_end;
        
    
        % LINE 12 deteting transmiiter activity
         y =  A(:,m);  
         y = zscore(y); 
    
        [W,~] = ADMM_spare_ortho_dic_encode(y',Phi',lambda,rho);
        W= W';
        W_c{i} = W;
    
        PhiW = Phi*W;
    
        % time step indeification
        options = statset('MaxIter',1000);
        NN = 2;
        gmm = fitgmdist(PhiW, NN, 'Options',options, ...
            'CovarianceType','diagonal', 'RegularizationValue',0.000001);
        cluster_means = gmm.mu;
        [~, mean_idx] = max(cluster_means);
        cluster_labels = cluster(gmm, PhiW);
    
        % building hole   
        ident_t=find(cluster_labels == mean_idx);
    
        ident_subs = zeros(length(ident_t)*length(ident_f),2);
        offset=1;
    
        for j=1:length(ident_t)
            ident_subs(offset:offset+length(ident_f)-1,1)=ident_t(j);
            ident_subs(offset:offset+length(ident_f)-1,2)=ident_f;
            offset  =offset+length(ident_f);
        end
    
        trans_ind = sub2ind([t,f] ,ident_subs(:,1),ident_subs(:,2));
    
        % stopping criterion
        det_box = X(trans_ind); 
        mean_det = (mean(det_box(:)));
        p_value = 1-normcdf(mean_det, pd.mu, pd.sigma);
    
        if (~isnan(p_value) && p_value <= threshold)
            pred{i} = zeros(t, f);
    
            Omega(trans_ind) = 0;
        
            pred{i}(trans_ind) = 1;
    
        else
            break;
    
        end
    
        i = i + 1;
    
    end


end


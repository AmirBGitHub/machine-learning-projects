function AdaBoost_Params() 

import = importdata('Digit Data.xlsx');
TrainingDataRead = import.data;

for c = 1:10
    for s = 1:10
        p = 1;
        if sum(isnan(TrainingDataRead(:,30*(c-1)+3*s-2))) == 0
            TrainingData{c,s} = TrainingDataRead(:,30*(c-1)+3*s-2:30*(c-1)+3*s);
        elseif sum(isnan(TrainingDataRead(:,30*(c-1)+3*s-2))) ~= 0
            while  isnan(TrainingDataRead(p,30*(c-1)+3*s-2)) == 0
                p = p + 1;
            end
            TrainingData{c,s} = TrainingDataRead(1:p-1,30*(c-1)+3*s-2:30*(c-1)+3*s);
        end    
        
    end
end

for c = 1:10
    for s = 1:10
        TrainingData_Features(10*(c-1) + s,:) = Stat_Geom_Features(TrainingData{c,s});
    end
end

for C1 = 1:10
    for C2 = C1+1:10
        % Training Data  
        % with 51 features for each sample 
        Cls1 = TrainingData_Features(10*(C1-1)+1:10*C1,:);
        Cls2 = TrainingData_Features(10*(C2-1)+1:10*C2,:);

        Features = [Cls1;Cls2];
        Class(1:10) = -1;
        Class(11:20) = 1;

        % Use Adaboost to make a classifier
        [ClassEstimate,model] = AdaBoost('train',Features,Class,20);
        Model{C1,C2} = model;
    end
end

save('TrainingParams.mat','Model')

%% Functions

    function [Class_Estimate,Params] = AdaBoost(Mode,Features,Class_or_model,Iteration)

        switch(Mode)
            case 'train' 

                Class = Class_or_model(:);    
                Params = struct;

                D = ones(length(Class),1)/length(Class);       
                Class_Estimate_Sum = zeros(size(Class));       
                boundary = [min(Features,[],1) max(Features,[],1)];        

                for t=1:Iteration
                    [Class_Estimate, err, h] = Weighted_Weak_Learners(Features,Class,D);      

                    alpha = 1/2 * log((1-err)/max(err,eps));   

                    Params(t).alpha = alpha;                   
                    Params(t).dimension = h.dimension;
                    Params(t).threshold = h.threshold;
                    Params(t).direction = h.direction;
                    Params(t).boundary = boundary;

                    D = (D.* exp(-Params(t).alpha.*Class.*Class_Estimate))./sum(D);        

                    Class_Estimate_Sum = Class_Estimate_Sum + Class_Estimate*Params(t).alpha;        
                    Class_Estimate = sign(Class_Estimate_Sum);
                    Params(t).error = sum(Class_Estimate ~= Class)/length(Class);
                    if(Params(t).error == 0) 
                        break; 
                    end
                end

            case 'apply'         

                Params = Class_or_model;    

                if(length(Params)>1);               
                    minb = Params(1).boundary(1:end/2);
                    maxb = Params(1).boundary(end/2+1:end);
                    Features = bsxfun(@min,Features,maxb);
                    Features = bsxfun(@max,Features,minb);
                end

                Class_Estimate_Sum = zeros(size(Features,1),1);
                for t=1:length(Params);
                    Class_Estimate_Sum = Class_Estimate_Sum + Params(t).alpha*Mu_Treshold_Check(Params(t), Features);
                end
                Class_Estimate = sign(Class_Estimate_Sum);

            otherwise
                error('Unknown Mode');
        end
    end


    function y = Mu_Treshold_Check(h, x)
        % Separates the data into two subsets based on the Mu and Feature values
        if(h.direction == 1)
            y =  double(x(:,h.dimension) >= h.threshold);
        else
            y =  double(x(:,h.dimension) < h.threshold);
        end
        y(y==0) = -1;

        end

    function [Class_Estimate,err,h] = Weighted_Weak_Learners(Features,Class,Weight)

        % Number of treshold steps
        ntre=2e5;

        % Split the data in two classes 1 and -1
        c1 = Features(Class<0,:); 
        w1 = Weight(Class<0);
        c2 = Features(Class>0,:); 
        w2 = Weight(Class>0);

        % Calculate the min and max for every dimensions
        mind = min(Features,[],1)-1e-10; 
        maxd = max(Features,[],1)+1e-10;

        % Make a weighted histogram of the two classes
        p2c = ceil((bsxfun(@rdivide,bsxfun(@minus,c2,mind),(maxd-mind)))*(ntre-1)+1+1e-9);   
        p2c(p2c>ntre) = ntre;
        p1f = floor((bsxfun(@rdivide,bsxfun(@minus,c1,mind),(maxd-mind)))*(ntre-1)+1-1e-9);  
        p1f(p1f<1) = 1;
        ndims = size(Features,2);
        i1 = repmat(1:ndims,size(p1f,1),1);  
        i2 = repmat(1:ndims,size(p2c,1),1);
        h1f = accumarray([p1f(:) i1(:)],repmat(w1(:),ndims,1),[ntre ndims],[],0);
        h2c = accumarray([p2c(:) i2(:)],repmat(w2(:),ndims,1),[ntre ndims],[],0);

        % This function calculates the error 
        h2ic = cumsum(h2c,1);
        h1rf = cumsum(h1f(end:-1:1,:),1); 
        h1rf = h1rf(end:-1:1,:);
        e1a = h1rf + h2ic;
        e2a = sum(Weight)-e1a;

        % Threshold value and dimension with the minimum error
        [err1a,ind1a] = min(e1a,[],1);  
        dim1a = (1:ndims); 
        dir1a = ones(1,ndims);
        [err2a,ind2a] = min(e2a,[],1);  
        dim2a = (1:ndims); 
        dir2a = -ones(1,ndims);
        A = [err1a(:),dim1a(:),dir1a(:),ind1a(:); err2a(:),dim2a(:),dir2a(:),ind2a(:)];
        [err,i] = min(A(:,1)); 
        dim = A(i,2); 
        dir = A(i,3); 
        ind = A(i,4);
        thresholds = linspace(mind(dim),maxd(dim),ntre);
        thr = thresholds(ind);

        % New threshold
        h.dimension = dim; 
        h.threshold = thr; 
        h.direction = dir;
        Class_Estimate = Mu_Treshold_Check(h,Features);

    end

end
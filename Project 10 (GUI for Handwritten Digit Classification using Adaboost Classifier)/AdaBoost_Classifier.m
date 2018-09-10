function index = AdaBoost_Classifier(INK_Data,Model)

INK_Data_Features = Stat_Geom_Features(INK_Data);

count = 1; 
indx = zeros(45,1);
for C1 = 1:10
    for C2 = C1+1:10

        % Read test data
        INK_Data= INK_Data_Features;
     
        % Classify the testdata using the trained model
         TestClass = AdaBoost('apply',INK_Data,Model{C1,C2});

        
         if TestClass==-1
           indx(count,1) = C1;
         elseif TestClass==1
           indx(count,1) = C2;
         end

        count = count + 1;

    end
end  

index = mode(indx);  

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

        % Threshold steps
        tre = 2e5;

        % Split the data in two classes 1 and -1
        c1 = Features(Class<0,:); 
        w1 = Weight(Class<0);
        c2 = Features(Class>0,:); 
        w2 = Weight(Class>0);

        % Calculate the min and max for every dimensions
        mind = min(Features,[],1)-1e-10; 
        maxd = max(Features,[],1)+1e-10;

        % Make a weighted histogram of the two classes
        p2c = ceil((bsxfun(@rdivide,bsxfun(@minus,c2,mind),(maxd-mind)))*(tre-1)+1+1e-9);   
        p2c(p2c>tre) = tre;
        p1f = floor((bsxfun(@rdivide,bsxfun(@minus,c1,mind),(maxd-mind)))*(tre-1)+1-1e-9);  
        p1f(p1f<1) = 1;
        ndims = size(Features,2);
        i1 = repmat(1:ndims,size(p1f,1),1);  
        i2 = repmat(1:ndims,size(p2c,1),1);
        h1f = accumarray([p1f(:) i1(:)],repmat(w1(:),ndims,1),[tre ndims],[],0);
        h2c = accumarray([p2c(:) i2(:)],repmat(w2(:),ndims,1),[tre ndims],[],0);

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
        M = [err1a(:),dim1a(:),dir1a(:),ind1a(:); err2a(:),dim2a(:),dir2a(:),ind2a(:)];
        [err,i] = min(M(:,1)); 
        dim = M(i,2); 
        dir = M(i,3); 
        ind = M(i,4);
        thresholds = linspace(minr(dim),maxr(dim),tre);
        thr = thresholds(ind);

        % New threshold
        h.dimension = dim; 
        h.threshold = thr; 
        h.direction = dir;
        Class_Estimate = Mu_Treshold_Check(h,Features);

    end

end
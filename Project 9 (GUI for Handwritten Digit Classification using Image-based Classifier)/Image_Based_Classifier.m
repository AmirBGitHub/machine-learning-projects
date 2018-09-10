function [score,index] = Image_Based_Classifier(INK_Data,TrainingData)


%% polar Transformation

for C=1:10
    Template{C} = TrainingData{1,C};
    Template{C} = Polar_Transform(Template{C});
end
Data = Polar_Transform(INK_Data);


%% Image Conversion

for C=1:10
    Template_img{C} = Image_Convert(Template{C},48);
end
Data_img = Image_Convert(Data,48);

%% Scores Finding

Score = zeros(10,4);
for C=1:10
    [HK, MHD, Tanimoto, Yule] = Classifiers(Data_img,Template_img{C});
    Score(C,:) = [HK, MHD, Tanimoto, Yule];
end

%% Normalizing

m = min(Score);
M = max(Score);
for C=1:10
    Score(C,:) = (Score(C,:)-m)./(M-m);
end


%% Classifying

[score,index] = Recognize(Score);

%% Functions

    function newPoints = Polar_Transform(p)
        
        l = sqrt(diff(p(:,1)).^2 + diff(p(:,2)).^2);
        
        xc = (sum(l'*p(1:end-1,1)) / sum(l) + sum(l'*p(2:end,1)) / sum(l))/2;   
        yc = (sum(l'*p(1:end-1,2)) / sum(l) + sum(l'*p(2:end,2)) / sum(l))/2;   
      
        r = sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2);
        
        p = (2*pi/max(r))*p;
        l = (2*pi/max(r))*l;    
        xc = (2*pi/max(r))*xc;
        yc = (2*pi/max(r))*yc;
        r = sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2);
        
        theta = atan2(p(:,2)-yc , p(:,1)-xc);
        newPoints = [theta,r];
    end



    function imag = Image_Convert(p,size)
        X = p(:,1);
        Y = p(:,2);
        H = max(Y)-min(Y);
        L = max(X)-min(X);
        d = abs(H-L)/2;
        
        if H>L
            X2 = (X+d/2)/(L+d)*47;
            Y2 = Y/H*47;
        else
            Y2 = (Y+d/2)/(H+d)*47;
            X2 = X/L*47;
        end
        
        X2 = X2-min(X2)+1;
        Y2 = Y2-min(Y2)+1;
        imag = zeros(size);
        for i=1:length(X2)
            imag(int8(X2(i)),int8(Y2(i))) = 1;
        end
    end
        


    function [HK,MHD,Tanimoto,Yule] = Classifiers(A,B)        
        
        % Pixel Locator
        n = 1;
        xy_a = [];
        for i = 1:size(A,1)
            for j = 1:size(A,2)
                if A(i,j)==1
                    xy_a(n,:)=[i,j];
                    n = n + 1;
                end
            end
        end
        
        n = 1;
        xy_b = [];
        for i = 1:size(B,1)
            for j = 1:size(B,2)
                if B(i,j)==1
                    xy_b(n,:)=[i,j];
                    n = n + 1;
                end
            end
        end
        
        
        % Classifier 1: HK
        for i = 1:length(xy_a)
            h_AB(i) = min(sqrt((xy_b(:,1) - xy_a(i,1)).^2 + (xy_b(:,2) - xy_a(i,2)).^2));
        end
        h_AB = sort(h_AB);
        HK_a = h_AB(round(length(xy_a)*0.94));
        
        for i = 1:length(xy_b)
            h_BA(i) = min(sqrt((xy_a(:,1) - xy_b(i,1)).^2 + (xy_a(:,2) - xy_b(i,2)).^2));
        end
        h_BA = sort(h_BA);
        HK_b = h_BA(round(length(xy_b)*0.94));
        
        HK = max(HK_a,HK_b);
        
        % Classifier 2: MHD
        for i = 1:length(xy_a)
            hmod_AB = (xy_a(i,2)^0.1)*sqrt((xy_b(:,1) - xy_a(i,1)).^2 + (xy_b(:,2) - xy_a(i,2)).^2);
        end
        MHD_a = mean(hmod_AB);
        
        for i = 1:length(xy_b)
            hmod_BA = (xy_b(i,2)^0.1)*sqrt((xy_a(:,1) - xy_b(i,1)).^2 + (xy_a(:,2) - xy_b(i,2)).^2);   
        end
        MHD_b = mean(hmod_BA);
        
        MHD = max(MHD_a,MHD_b);
        
        % Classifier 3 & 4
        na = sum(sum(A==1));
        nb = sum(sum(B==1));
        n00 = sum(sum((1-A)+(1-B)==2));

        % Overlapping
        nab = 0;
        dist = pdist2(xy_a,xy_b);
        for i=1:length(xy_a)
            if min(dist(i,:))<4.5
                nab = nab+1;
            end
        end
        
        T_AB = nab/(na + nb - nab);
        TC_AB = n00/(na + nb - 2*nab + n00);
        
        a = 0.75 - 0.25*((na + nb)/(2*48*48));
        
        % Tonimoto
        Tanimoto = a*T_AB + (1-a)*TC_AB;
        
        % Yule
        Yule = (nab*n00 - (na - nab)*(nb - nab))/(nab*n00 + (na - nab)*(nb - nab));
        
    end



    function [score,index] = Recognize(Score)
        [V(1),I(1)] = min(Score(:,1));
        [V(2),I(2)] = min(Score(:,2));
        [V(3),I(3)] = max(Score(:,3));
        [V(4),I(4)] = max(Score(:,4));
        index = mode(I);
        score = I;
    end

end
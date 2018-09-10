function Features = Stat_Geom_Features(Points)
    
    X = zeros(51,1);
    n = length(Points);

    %% Symbol strokes
    X(1,:) = sum(Points(:,3)); 
    EndPntIndx = [find(Points(:,3)); n-1];
    
    %% Cusp Features
    d = diff(Points(:,1:2));
    CuspNum = 0;
    CuspIndex = [];
    for i = 1:n-2
        alpha = acosd(d(i,:)*d(i+1,:)'/(norm(d(i,:))*norm(d(i+1,:))));
        if alpha > 60
            CuspIndex = [CuspIndex; i];
            CuspNum = CuspNum + 1;
        end
    end
    X(2,:) = CuspNum;
    if CuspNum ~= 0
        X(3,:) = max(sqrt(sum((Points(CuspIndex,1:2)-repelem(Points(1,1:2),CuspNum,1)).^2,2)));
        X(4,:) = min(sqrt(sum((Points(CuspIndex,1:2)-repelem(Points(1,1:2),CuspNum,1)).^2,2)));
        X(5,:) = max(sqrt(sum((Points(CuspIndex,1:2)-repelem(Points(end,1:2),CuspNum,1)).^2,2)));
        X(6,:) = min(sqrt(sum((Points(CuspIndex,1:2)-repelem(Points(end,1:2),CuspNum,1)).^2,2)));   
    else
        X(3,:) = 0;
        X(4,:) = 0;
        X(5,:) = 0;
        X(6,:) = 0;
    end
    
    %% Aspect ratio
    Bounding_Box = [max(Points(:,1)) - min(Points(:,1)); max(Points(:,2)) - min(Points(:,2))];
    X(7,:) = Bounding_Box(1)/Bounding_Box(2);
    
    %% intersection features
    IntersectNum = 0;
    IntersectPnts = [];
    syms x y
    
    for i = 1:n-3
        for j = i+2:n-1
            if sign(det([Points(i+1,1) Points(i+1,2) 1; Points(j,1) Points(j,2) 1; Points(j+1,1) Points(j+1,2) 1])*det([Points(i,1) Points(i,2) 1; Points(j,1) Points(j,2) 1; Points(j+1,1) Points(j+1,2) 1])) < 0 && sign(det([Points(i,1) Points(i,2) 1; Points(i+1,1) Points(i+1,2) 1; Points(j,1) Points(j,2) 1])*det([Points(i,1) Points(i,2) 1; Points(i+1,1) Points(i+1,2) 1; Points(j+1,1) Points(j+1,2) 1]))<0 
                IntersectNum = IntersectNum + 1;
%                 [soly, solx] = solve([ y - ((Points(i+1,2)-Points(i,2))/(Points(i+1,1)-Points(i,1)))*x + Points(i,1)*((Points(i+1,2)-Points(i,2))/(Points(i+1,1)-Points(i,1))) - Points(i,2) == 0,  y - ((Points(j+1,2)-Points(j,2))/(Points(j+1,1)-Points(j,1)))*x +  Points(j,1)*((Points(j+1,2)-Points(j,2))/(Points(j+1,1)-Points(j,1))) - Points(j,2) == 0], [y, x]);
%                 IntersectPnts = [IntersectPnts; soly, solx]     
            end
        end 
    end
    X(8,:) = IntersectNum;
%     if IntersectNum ~= 0 || isempty(IntersectPnts) == 0 
%         X(9,:) = max(sqrt(sum((IntersectPnts-repelem(Points(1,1:2),IntersectNum,1)).^2,2)));
%         X(10,:) = min(sqrt(sum((IntersectPnts-repelem(Points(1,1:2),IntersectNum,1)).^2,2)));
%         X(11,:) = max(sqrt(sum((IntersectPnts-repelem(Points(end,1:2),IntersectNum,1)).^2,2)));
%         X(12,:) = min(sqrt(sum((IntersectPnts-repelem(Points(end,1:2),IntersectNum,1)).^2,2)));
%     else
        X(9,:) = 0;
        X(10,:) = 0;
        X(11,:) = 0;
        X(12,:) = 0;
%    end
    
    %% 2D point histogram
    xl = min(Points(:,1));
    yl = min(Points(:,2));
    dx = Bounding_Box(1)/3;
    dy = Bounding_Box(2)/3;
    
    PntHist = zeros(3);
    for i = 0:2
        for j = 0:2
            PntHist(i+1,j+1) = numel(find(Points(:,1) > xl+i*dx & Points(:,1) < xl+(i+1)*dx & Points(:,2) > yl+j*dy & Points(:,2) < yl+(j+1)*dy))/numel(Points(:,1));
        end
    end
    
    for i = 1:3
        for j = 1:3
            X(12+3*(i-1)+j,:) = PntHist(i,j);
        end
    end
    
    %% Angle histogram
    X_p = zeros(X(1,:),8);
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        J = length(StrokePnts)-1;
        AngHist = zeros(1,8);
        for j = 1:J
            v = diff(StrokePnts);
            theta = acosd([1 0]*(v/norm(v))')';  
        end
        for k = 0:7
            AngHist(:,k+1) = numel(find(theta > k*45 & theta <= (k+1)*45))/numel(theta);
        end
        X_p(i,:) = AngHist;
    end
    X(22:29,:) = mean(X_p,2);
    
    %% First and last distance
    dist = zeros(X(1,:),1);
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        dist(i,:) = norm(StrokePnts(end,:) - StrokePnts(1,:));
    end
    X(30,:) = mean(dist);
    
    %% Arc length
    X_p = zeros(X(1,:),1);
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        X_p(i,:) = sum(sqrt(sum(diff(StrokePnts).^2,2)));
    end
    X(31,:) = sum(X_p);
    
    %% Fit line feature
    err = zeros(X(1,:),1);
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        p = polyfit(StrokePnts(:,1),StrokePnts(:,2),1);
        err(i,:) = sum((StrokePnts(:,2)-polyval(p,StrokePnts(:,1))).^2);
    end
    X(32,:) = sum(err);
    
    
    %% Dominant point feature
    a_max = zeros(X(1,:),1);
    a_dev = zeros(X(1,:),1);
    a_cnt = zeros(X(1,:),1);
    a_zc = zeros(X(1,:),1);
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        m = length(StrokePnts);
        d = diff(StrokePnts(:,1:2));
        for i = 4:m-4
            %betha(i,:) = mean( acosd(d(i,:)*d(i+1,:)'/(norm(d(i,:))*norm(d(i+1,:)))) + acosd((d(i-1,:)+d(i,:))*(d(i+1,:)+d(i+2,:))'/(norm(d(i-1,:)+d(i,:))*norm(d(i+1,:)+d(i+2,:)))) +  acosd((d(i-2,:)+d(i-1,:)+d(i,:))*(d(i+1,:)+d(i+2,:)+d(i+3,:))'/(norm(d(i-2,:)+d(i-1,:)+d(i,:))*norm(d(i+1,:)+d(i+2,:)+d(i+3,:)))) + acosd((d(i-3,:)+d(i-2,:)+d(i-1,:)+d(i,:))*(d(i+1,:)+d(i+2,:)+d(i+3,:)+d(i+4,:))'/(norm((d(i-3,:)+d(i-2,:)+d(i-1,:)+d(i,:))*norm(d(i+1,:)+d(i+2,:)+d(i+3,:)+d(i+4,:))))) ); 
            betha(i,:) = acosd((d(i-3,:)+d(i-2,:)+d(i-1,:)+d(i,:))*(d(i+1,:)+d(i+2,:)+d(i+3,:))'/(norm((d(i-3,:)+d(i-2,:)+d(i-1,:)+d(i,:))*norm(d(i+1,:)+d(i+2,:)+d(i+3,:)))));
        end
        [~, ExtrmLocs] = findpeaks(betha);
        MidPnts = round(mean([1 ExtrmLocs(1:end)';ExtrmLocs(1:end)' m]));
        DominantPnts = [StrokePnts(1,:); StrokePnts(ExtrmLocs,:); StrokePnts(MidPnts,:); StrokePnts(end,:)];
        mdp = length(DominantPnts);
        v = zeros(mdp-1,2);
        for j = 1:mdp-1
            v = diff(DominantPnts);
            gamma = acosd([1 0]*(v/norm(v))')'; 
        end
        a_max(i,:) = max(gamma);
        Phi = diff(gamma);
        a_dev(i,:) = (1/(mdp-2))*sum(Phi);
        a_cnt(i,:) = nnz(Phi<3)/(mdp-2);
        a_zc(i,:) = nnz(diff(sign(Phi(find(Phi)))));
    end
    X(33,:) = mean(a_max);
    X(34,:) = mean(a_dev);
    X(35,:) = mean(a_cnt);
    X(36,:) = mean(a_zc);
    
    %% Stroke area
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        m = length(StrokePnts);
        Triangle_area = zeros(m,1);
        for j = 1:m-2
            u = (StrokePnts(j+1,:) - StrokePnts(1,:));
            v = (StrokePnts(j+2,:) - StrokePnts(1,:));
            Triangle_area(j,:) = 0.5*det([u; v])*sign(det([u; v]));
        end
        S_area = sum(Triangle_area);
    end
    X(37,:) = mean(S_area);
    
    %% Side ratios
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        m = length(StrokePnts);
        xl = min(StrokePnts(:,1));
        Bounding_Box = [max(StrokePnts(:,1)) - min(StrokePnts(:,1)); max(StrokePnts(:,2)) - min(StrokePnts(:,2))];
        X_p_1(i,:) = (StrokePnts(1,1) - xl)/Bounding_Box(1);
        X_p_2(i,:) = (StrokePnts(end,1) - xl)/Bounding_Box(1);
    end
    X(38,:) = mean(X_p_1);
    X(39,:) = mean(X_p_2);
    
    %% Top and bottom ratios
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        yl = max(StrokePnts(:,2));
        Bounding_Box = [max(StrokePnts(:,1)) - min(StrokePnts(:,1)); max(StrokePnts(:,2)) - min(StrokePnts(:,2))];
        X_p_1(i,:) = (StrokePnts(1,2) - yl)/Bounding_Box(2);
        X_p_2(i,:) = (StrokePnts(end,2) - yl)/Bounding_Box(2);
    end
    X(40,:) = mean(X_p_1);
    X(41,:) = mean(X_p_2);
    
    %% Min and max features
    for i = 1:X(1,:)
        StrokePnts = Points(EndPntIndx(i):EndPntIndx(i+1)-1,1:2);
        Start_dir_x = zeros(X(1,:),1);
        End_dir_x = zeros(X(1,:),1);
        Start_dir_y = zeros(X(1,:),1);
        End_dir_y = zeros(X(1,:),1);
        
        if (StrokePnts(2,1)-StrokePnts(1,1))<=0 && (StrokePnts(end,1)-StrokePnts(end-1,1))>=0
            Start_dir_x(i,:) = -1;
            End_dir_x(i,:) = 1;
            [~, min_locs] = findpeaks(-StrokePnts(:,1)); 
            local_minima_x = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,1)); 
            local_maxima_x = numel(max_locs);
            if local_minima_x == 0 || local_maxima_x == 0
                last_dis_x = 0;
            elseif local_minima_x ~= 0 || local_maxima_x ~= 0
                last_dis_x = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        elseif (StrokePnts(2,1)-StrokePnts(1,1))>=0 && (StrokePnts(end,1)-StrokePnts(end-1,1))>=0
            Start_dir_x(i,:) = 1;
            End_dir_x(i,:) = 1;
            [~, min_locs] = findpeaks(-StrokePnts(:,1)); 
            local_minima_x = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,1)); 
            local_maxima_x = numel(max_locs);
            if local_minima_x == 0 || local_maxima_x == 0
                last_dis_x = 0;
            elseif local_minima_x ~= 0 || local_maxima_x ~= 0
                last_dis_x = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        elseif (StrokePnts(2,1)-StrokePnts(1,1))<=0 && (StrokePnts(end,1)-StrokePnts(end-1,1))<=0
            Start_dir_x(i,:) = -1;
            End_dir_x(i,:) = -1;
            [~, min_locs] = findpeaks(-StrokePnts(:,1)); 
            local_minima_x = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,1)); 
            local_maxima_x = numel(max_locs);
            if local_minima_x == 0 || local_maxima_x == 0
                last_dis_x = 0;
            elseif local_minima_x ~= 0 || local_maxima_x ~= 0
                last_dis_x = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        elseif (StrokePnts(2,1)-StrokePnts(1,1))>=0 && (StrokePnts(end,1)-StrokePnts(end-1,1))<=0
            Start_dir_x(i,:) = 1;
            End_dir_x(i,:) = -1;
            [~, min_locs] = findpeaks(-StrokePnts(:,1)); 
            local_minima_x = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,1)); 
            local_maxima_x = numel(max_locs);
            if local_minima_x == 0 || local_maxima_x == 0
                last_dis_x = 0;
            elseif local_minima_x ~= 0 || local_maxima_x ~= 0
                last_dis_x = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        end 
        
        if (StrokePnts(2,2)-StrokePnts(1,2))<=0 && (StrokePnts(end,2)-StrokePnts(end-1,2))>=0
            Start_dir_y(i,:) = -1;
            End_dir_y(i,:) = 1;
            [~, min_locs] = findpeaks(-StrokePnts(:,2)); 
            local_minima_y = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,2)); 
            local_maxima_y = numel(max_locs);
            if local_minima_y == 0 || local_maxima_y == 0
                last_dis_y = 0;
            elseif local_minima_y ~= 0 || local_maxima_y ~= 0
                last_dis_y = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        elseif (StrokePnts(2,2)-StrokePnts(1,2))>=0 && (StrokePnts(end,2)-StrokePnts(end-1,2))>=0
            Start_dir_y(i,:) = 1;
            End_dir_y(i,:) = 1;
            [~, min_locs] = findpeaks(-StrokePnts(:,2)); 
            local_minima_y = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,2)); 
            local_maxima_y = numel(max_locs);
            if local_minima_y == 0 || local_maxima_y == 0
                last_dis_y = 0;
            elseif local_minima_y ~= 0 || local_maxima_y ~= 0
                last_dis_y = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        elseif (StrokePnts(2,2)-StrokePnts(1,2))<=0 && (StrokePnts(end,2)-StrokePnts(end-1,2))<=0
            Start_dir_y(i,:) = -1;
            End_dir_y(i,:) = -1;
            [~, min_locs] = findpeaks(-StrokePnts(:,2)); 
            local_minima_y = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,2)); 
            local_maxima_y = numel(max_locs);
            if local_minima_y == 0 || local_maxima_y == 0
                last_dis_y = 0;
            elseif local_minima_y ~= 0 || local_maxima_y ~= 0
                last_dis_y = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        elseif (StrokePnts(2,2)-StrokePnts(1,2))>=0 && (StrokePnts(end,2)-StrokePnts(end-1,2))<=0
            Start_dir_y(i,:) = 1;
            End_dir_y(i,:) = -1;
            [~, min_locs] = findpeaks(-StrokePnts(:,2)); 
            local_minima_y = numel(min_locs);
            [~, max_locs] = findpeaks(StrokePnts(:,2)); 
            local_maxima_y = numel(max_locs);
            if local_minima_y == 0 || local_maxima_y == 0
                last_dis_y = 0;
            elseif local_minima_y ~= 0 || local_maxima_y ~= 0
                last_dis_y = norm(StrokePnts(end,:)-StrokePnts(max([min_locs; max_locs]),:));
            end
        end
        X_p_1(i,:) = sign(Start_dir_x);
        X_p_2(i,:) = sign(End_dir_x);
        X_p_3(i,:) = round(local_minima_x);
        X_p_4(i,:) = round(local_maxima_x);
        X_p_5(i,:) = last_dis_x;
        X_p_6(i,:) = sign(Start_dir_y);
        X_p_7(i,:) = sign(End_dir_y);
        X_p_8(i,:) = round(local_minima_y);
        X_p_9(i,:) = round(local_maxima_y);
        X_p_10(i,:) = last_dis_y;
        
    end
    
    X(42,:) = mean(X_p_1);
    X(43,:) = mean(X_p_2);
    X(44,:) = mean(X_p_3);
    X(45,:) = mean(X_p_4);
    X(46,:) = mean(X_p_5);
    X(47,:) = mean(X_p_6);
    X(48,:) = mean(X_p_7);
    X(49,:) = mean(X_p_8);
    X(50,:) = mean(X_p_9);
    X(51,:) = mean(X_p_10);
    
    Features = X';
end


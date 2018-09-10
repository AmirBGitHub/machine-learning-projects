clear all
close all
clc

TrainingDataRead = xlsread('Digit Data'); 

for c = 1:10
    for e = 1:10
        p = 1;
        if sum(isnan(TrainingDataRead(:,30*(c-1)+3*e-2))) == 0
            TrainingData{c,e} = TrainingDataRead(:,30*(c-1)+3*e-2:30*(c-1)+3*e);
        elseif sum(isnan(TrainingDataRead(:,30*(c-1)+3*e-2))) ~= 0
            while  isnan(TrainingDataRead(p,30*(c-1)+3*e-2)) == 0
                p = p + 1;
            end
            TrainingData{c,e} = TrainingDataRead(1:p-1,30*(c-1)+3*e-2:30*(c-1)+3*e);
        end    
        
    end
end

for c = 1:10
    for e = 1:10
        x{c,e} = TrainingData{c,e}(:,1);
        y{c,e} = TrainingData{c,e}(:,2);
        t{c,e} = TrainingData{c,e}(:,3);
        P{c,e} = length(TrainingData{c,e}(:,1));
    end
end

f = zeros(10,10,13);

for c = 1:10
    for e = 1:10
        
        f(c,e,1) = (x{c,e}(3)-x{c,e}(1))/sqrt((x{c,e}(3)-x{c,e}(1))^2+(y{c,e}(3)-y{c,e}(1))^2);
        f(c,e,2) = (y{c,e}(3)-y{c,e}(1))/sqrt((x{c,e}(3)-x{c,e}(1))^2+(y{c,e}(3)-y{c,e}(1))^2);
        f(c,e,3) = sqrt((max(x{c,e}(:))-min(x{c,e}(:)))^2+(max(y{c,e}(:))-min(y{c,e}(:)))^2);
        f(c,e,4) = atan2((max(y{c,e}(:))-min(y{c,e}(:)))/(max(x{c,e}(:))-min(x{c,e}(:))),1);
        f(c,e,5) = sqrt((x{c,e}(P{c,e})-x{c,e}(1))^2+(y{c,e}(P{c,e})-y{c,e}(1))^2);
        f(c,e,6) = (x{c,e}(P{c,e})-x{c,e}(1))/f(c,e,5);
        f(c,e,7) = (y{c,e}(P{c,e})-y{c,e}(1))/f(c,e,5);
        f(c,e,8) = 0;
        f(c,e,9) = 0;
        f(c,e,10) = 0;
        f(c,e,11) = 0;
        for p = 1:P{c,e}-1
            dx(p) = x{c,e}(p+1)-x{c,e}(p);
            dy(p) = y{c,e}(p+1)-y{c,e}(p);
            dt(p) = t{c,e}(p+1)-t{c,e}(p);
            f(c,e,8) = f(c,e,8) + sqrt(dx(p)^2+dy(p)^2); 
        end
        for p = 2:P{c,e}-1
            f(c,e,9) = f(c,e,9) + atan2((dx(p)*dy(p-1)-dx(p-1)*dy(p))/(dx(p)*dx(p-1)+dy(p)*dy(p-1)),1);
            f(c,e,10) = f(c,e,10) + abs(atan2((dx(p)*dy(p-1)-dx(p-1)*dy(p))/(dx(p)*dx(p-1)+dy(p)*dy(p-1)),1));
            f(c,e,11) = f(c,e,11) + (atan2((dx(p)*dy(p-1)-dx(p-1)*dy(p))/(dx(p)*dx(p-1)+dy(p)*dy(p-1)),1))^2;
        end    

        f(c,e,12) = max((dx(1:P{c,e}-1).^2+dy(1:P{c,e}-1).^2)./dt(1:P{c,e}-1).^2);
        f(c,e,13) = t{c,e}(P{c,e})-t{c,e}(1);
    end
end

for c = 1:10
    for i = 1:13   
        f_mean(c,i) = mean(f(c,:,i));
    end
end

for c = 1:10
    for i = 1:13
        for j = 1:13
            Cov(c,i,j) = [f(c,:,i) - f_mean(c,i)]*[f(c,:,j) - f_mean(c,j)]';
        end
    end
end

for i = 1:13
    for j = 1:13
        Avg_Cov(i,j) = sum(Cov(:,i,j))/90;
    end
end

for c = 1:10
    for j = 1:13
        Avg_Cov_inv = inv(Avg_Cov);
        w(c,j) = sum(Avg_Cov_inv(:,j)'.*f_mean(c,:));
    end
end

for c = 1:10
    w0(c,:) = -0.5*sum(w(c,:).*f_mean(c,:));
end

save('w.mat', 'w');
save('w0.mat', 'w0');




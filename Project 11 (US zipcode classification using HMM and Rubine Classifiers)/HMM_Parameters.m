close all
clear 
clc

zipcode_data = importdata('zipcode_data.mat');
zipcode = zipcode_data.zipcode;
state = zipcode_data.state;
W = importdata('W.mat');
w = W.w;
w0 = W.w0;
state1 = importdata('data_1.mat');
feature1 = state1.Fi;
Score(:,1) = w'*feature1 + w0';
state2 = importdata('data_2.mat');
feature2 = state2.Fi;
Score(:,2) = w'*feature2 + w0';
state3 = importdata('data_3.mat');
feature3 = state3.Fi;
Score(:,3) = w'*feature3 + w0';
state4 = importdata('data_4.mat');
feature4 = state4.Fi;
Score(:,4) = w'*feature4 + w0';
state5 = importdata('data_5.mat');
feature5 = state5.Fi;
Score(:,5) = w'*feature5 + w0';

%% NY

DgtTrnsCnt = zeros(10,10,4); 
DgtTrnsCnt_NY_sub = zeros(10,4);
strdigits = num2str(zipcode(1:length(state(state(:,1)=='N' & state(:,2)=='Y')),1));
d1 = str2num(strdigits(:,1));
d2 = str2num(strdigits(:,2));
d3 = str2num(strdigits(:,3));
d4 = str2num(strdigits(:,4));
d5 = str2num(strdigits(:,5));
zipcode_digits_NY = [d1, d2, d3, d4, d5];

% PIi
pi_NY(:,1) = [1 0 0 0 0 0 0 0 0 0]';

for k = 1:4
    for i = 1:10
        for j = 1:10
            DgtTrnsCnt_NY(i,j,k) = length(zipcode_digits_NY(zipcode_digits_NY(:,k) == i-1 & zipcode_digits_NY(:,k+1) == j-1,:));        
        end
    end
end

for k = 1:4
    for i = 1:10
        DgtTrnsCnt_NY_sub(i,k) = length(zipcode_digits_NY(zipcode_digits_NY(:,k) == i-1,:));        
    end
end


% aij
aij_NY = sum(DgtTrnsCnt_NY(:,:,:),3)./repmat(sum(DgtTrnsCnt_NY_sub(:,:),2),1,10);

% bj
for i = 1:10
    bj_NY(i,:) = Score(i,:).*(1./sum(Score));
end

save('NY_HMM_Parameters.mat', 'pi_NY', 'aij_NY', 'bj_NY');

%% US

DgtTrnsCnt = zeros(10,10,4); 
DgtTrnsCnt_US_sub = zeros(10,4);
strdigits = num2str(zipcode(:,1));
d1 = str2num(strdigits(:,1));
d2 = str2num(strdigits(:,2));
d3 = str2num(strdigits(:,3));
d4 = str2num(strdigits(:,4));
d5 = str2num(strdigits(:,5));
zipcode_digits_US = [d1, d2, d3, d4, d5];

% PIi
p1 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 0,:))/length(zipcode_digits_US); 
p2 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 1,:))/length(zipcode_digits_US); 
p3 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 2,:))/length(zipcode_digits_US); 
p4 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 3,:))/length(zipcode_digits_US); 
p5 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 4,:))/length(zipcode_digits_US); 
p6 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 5,:))/length(zipcode_digits_US); 
p7 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 6,:))/length(zipcode_digits_US); 
p8 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 7,:))/length(zipcode_digits_US); 
p9 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 8,:))/length(zipcode_digits_US); 
p0 = length(zipcode_digits_US(zipcode_digits_US(:,1) == 9,:))/length(zipcode_digits_US); 

pi_US(:,1) = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9]';

for k = 1:4
    for i = 1:10
        for j = 1:10
            DgtTrnsCnt_US(i,j,k) = length(zipcode_digits_US(zipcode_digits_US(:,k) == i-1 & zipcode_digits_US(:,k+1) == j-1,:));        
        end
    end
end

for k = 1:4
    for i = 1:10
        DgtTrnsCnt_US_sub(i,k) = length(zipcode_digits_US(zipcode_digits_US(:,k) == i-1,:));        
    end
end


% aij
aij_US = sum(DgtTrnsCnt_US(:,:,:),3)./repmat(sum(DgtTrnsCnt_US_sub(:,:),2),1,10);

% bj
for i = 1:10
    bj_US(i,:) = Score(i,:).*(1./sum(Score));
end

save('US_HMM_Parameters.mat', 'pi_US', 'aij_US', 'bj_US');






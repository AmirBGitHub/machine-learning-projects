close all
clear
clc

US_ProbData = load('US_HMM_Parameters.mat');
pi_US = US_ProbData.pi_US;
aij_US = US_ProbData.aij_US;
bj_US = US_ProbData.bj_US;

Classes = [0 1 2 3 4 5 6 7 8 9]';
Output = [1 4 2 6 3];

delta = zeros(5,10);
path = zeros(1,5);

delta(1,:) = pi_US.*bj_US(:,1);
[Pi, maxima(1)] = max(delta(1,:));


for i = 2:5
    for j = 1:10
        [delta(i,j), maxim(i,j)] = max(delta(i-1,:)'.*aij_US(:,j).*bj_US(j,i));
    end  
end

[V, state_class(5)] = max(delta(5,:));
for i = 4:-1:1
    state_class(i) = maxim(i+1,state_class(i+1));
end

Zipcode_Classified = Output(state_class);
save('US_HMM_Classified.mat', 'Zipcode_Classified');



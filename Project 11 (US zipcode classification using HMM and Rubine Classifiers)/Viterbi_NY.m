close all
clear
clc

NY_ProbData = load('NY_HMM_Parameters.mat');
pi_NY = NY_ProbData.pi_NY;
aij_NY = NY_ProbData.aij_NY;
bj_NY = NY_ProbData.bj_NY;

Classes = [0 1 2 3 4 5 6 7 8 9]';
Output = [1 4 2 6 3];

delta = zeros(5,10);
path = zeros(1,5);

delta(1,:) = pi_NY.*bj_NY(:,1);
[Pi, maxima(1)] = max(delta(1,:));


for i = 2:5
    for j = 1:10
        [delta(i,j), maxim(i,j)] = max(delta(i-1,:)'.*aij_NY(:,j).*bj_NY(j,i));
    end  
end

[V, state_class(5)] = max(delta(5,:));
for i = 4:-1:1
    state_class(i) = maxim(i+1,state_class(i+1));
end

Zipcode_Classified = Output(state_class);
save('NY_HMM_Classified.mat', 'Zipcode_Classified');


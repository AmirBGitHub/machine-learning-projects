close all
clear
clc

W = importdata('W.mat');
w = W.w;
w0 = W.w0;
state1 = importdata('data_1.mat');
feature1 = state1.Fi;
[v1 c1] = max(w'*feature1 + w0');
class1 = mod(c1,10);

state2 = importdata('data_2.mat');
feature2 = state2.Fi;
[v2 c2] = max(w'*feature2 + w0');
class2 = mod(c2,10);

state3 = importdata('data_3.mat');
feature3 = state3.Fi;
[v3 c3] = max(w'*feature3 + w0');
class3 = mod(c3,10);

state4 = importdata('data_4.mat');
feature4 = state4.Fi;
[v4 c4] = max(w'*feature4 + w0');
class4 = mod(c4,10);

state5 = importdata('data_5.mat');
feature5 = state5.Fi;
[v5 c5] = max(w'*feature5 + w0');
class5 = mod(c5,10);

Zipcode_Classified = [class1, class2, class3, class4, class5];

save('Rubine_Classified.mat', 'Zipcode_Classified');


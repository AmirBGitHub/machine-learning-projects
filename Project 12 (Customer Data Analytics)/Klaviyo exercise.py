# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 10:39:56 2018

@author: Amir Baghdadi
Email: amirbagh@buffalo.edu / amir.baqdadi@gmail.com'
"""

import csv
import pandas as pd
import numpy as np
from datetime import date
from collections import Counter
import matplotlib.pyplot as plt
from scipy import stats
from IPython.display import display

print(' ')
print('By: Amir Baghdadi')
print('Email: amirbagh@buffalo.edu / amir.baqdadi@gmail.com')
print(' ')

data_read = []
with open('data_science_screening_exercise_orders.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        data_read.append(row) 

data = pd.DataFrame(data_read)
# sorting the data based on customer ID keeping the old indices
data_sort_oldidx = data.sort_values(by=[0])
# sorting the data based on customer ID updating the indices                 
data_sort_newidx = data_sort_oldidx.reset_index(drop=True)  

# making a dataframe of IDs along with the order counts 
IDs_list = data_sort_newidx[0].tolist()
IDs_dataframe = pd.DataFrame.from_dict({x:IDs_list.count(x) for x in IDs_list}, orient='index')

###############################################################################  
# assembling a dataframe with one row per customer and the following columns:
#    * customer_id
#    * gender
#    * most_recent_order_date
#    * order_count (number of orders placed by this customer)
customer_list = []
for i in range(len(IDs_dataframe)):   
    customer_list.append([IDs_dataframe.index[i], \
                          data_sort_newidx[1][sum(IDs_dataframe[0][0:i])], \
                          data_sort_oldidx[2][max(data_sort_oldidx.index[range(sum(IDs_dataframe[0][0:i]),sum(IDs_dataframe[0][0:i+1]))])], \
                          IDs_dataframe[0][i]])

customer_dataframe = pd.DataFrame(customer_list).rename(columns = {0:'customer_id', \
                                                                   1:'gender', \
                                                                   2:'most_recent_order_date', \
                                                                   3:'order_count'})
print("A)")                                                                  
print("Customer Orders Dataframe:")
display(customer_dataframe)
customer_dataframe.to_csv("customer orders dataframe.csv", sep=',')
###############################################################################
print("B)")      
# plotting data per week
# starting the week from 1/1/2017 and considering each consecutive 7 days as a week    
origin = date(2017, 1, 1)
weekNum = []       
for i in range(len(data_read)-1):        
    dateTime = data_read[i+1][2].split(' ')[0].split('-') 
    order_date = date(int(dateTime[0]),int(dateTime[1]),int(dateTime[2]))
    weekNum.append(int((order_date - origin).days / 7 + 1))
    
weeklyStat = Counter(weekNum)
plt.plot(weeklyStat.keys(), weeklyStat.values(), 'b--')
plt.xlabel('Week Number')
plt.ylabel('Order Counts')
plt.title('Count of Orders per Week in 2017')
plt.savefig('data per week plot')
plt.show()        

###############################################################################
print("C)")          
# running pairwise t-test for comparing the mean order value for gender 0 and 1    
gender_0_values = pd.to_numeric(data.loc[data[1] == '0'][3])
mean_order_value_gender_0 = np.mean(gender_0_values)
gender_1_values = pd.to_numeric(data.loc[data[1] == '1'][3])
mean_order_value_gender_1 = np.mean(gender_1_values)
t_test = stats.ttest_ind(gender_0_values,gender_1_values)
print(t_test)
p_value = round(t_test[1],3)
print('p-value = ', p_value)
# p-value < 0.05 so the difference between the means is significant however very marginal
print('p-value < 0.05 so the difference between the means is significant however very marginal') 



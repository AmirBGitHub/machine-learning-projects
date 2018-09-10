# machine-learning-projects
My machine learning and deep learning projects

# Guideline

Projects 1-4: After cloning/downloading the repository and extracting the ZIP file and executing command python main.py in the first level directory, it will be able to generate all the results and plots used in the PDF report and print them out in a clear manner.

Project 5: After cloning/downloading the repository and running the 'alphaBuildFeatures.m' file, it will be able to generate the results in two separate '.mat' files. The classification codes and results for surgeon robotic skill classification are available in the 'classification results' folder.

Project 6-10: After cloning/downloading the repository right click on 'INK.fig' from MATLAB and 'Open in GUIDE'. After running the GUI  your handwritten curve will be segmented or the digit will be classified.

Project 11: After cloning/downloading the repository run the 'Rubine.m', 'Viterbi_NY.m', or, 'Viterbi_US.m' to get the different zipcode classification results.

Project 12: After cloning/downloading the repository run the 'Klaviyo exercise.py' file to get the statistical analysis results for customer data.

# Project 1

The purpose of this project was to calculate some statistics and probability concepts using the data of ranking for US universities in CS with four criteria of CS ranking based on US news, the value of Research Overhead, Admin Base Pay, and Out-of-state Tuition. 

In addition, the probability of the variables have been determined and, by considering the correlation between each variable and constructing a Bayesian network, this probability has been improved by obtaining a higher probability for the variables in terms of log-likelihood.

# Project 2

The purpose of this project is to use machine learning approach to solve a LeToR problem. The problem is formulated as a linear regression where the input vector x is mapped to a real-valued scalar target y(x, w). The tasks included training two LeToR and synthetic datasets using closedform and stochastic gradient descent (SGD) solutions. 

The data consisted of pairs of input values x with different number of features, i.e, 10 for synthetic and 46 for LeToR data, and target value t. The input values are real-valued vectors (features derived from a query-document pair) and the target values are scalars (relevance labels) that take one of three values 0, 1, 2, where the larger relevance label corresponds with the better match between query and document. Despite of having discrete targets, linear regression was used to obtain real values that avoids restricting to only three possible values of targets which is more useful in ranking.

# Project 3

The purpose of this project is to implement a multiclass LR and a single hidden layer ANN classification algorithms on MNIST dataset and evaluate the performance on both MNIST and USPS handwritten digit data to identify them among 0, 1, 2, … , 9. The MNIST database is a large database of handwritten digits that is commonly used for training various image processing systems. 

The database contains 60,000 training images and 10,000 testing images. The images are centered in a 28x28 image by computing the center of mass of the pixels, and translating the image so as to position this point at the center of the 28x28 field. Each digit in the USPS dataset has 150 samples of “.png” file available for testing that need to be imported and converted into MNIST data format. The evaluation will be performed by calculating classification error rate in under the one-hot coding scheme.

# Project 4

The purpose of this project is to implement convolutional convolution neural network to determine whether the person in a portrait image is wearing glasses or not. The Celeb dataset was used which has more than 200k celebrity images in total. A publicly available convolutional neural network package from Tensorflow was used and trained on the Celeb images, hyperparameters were tuned and regularization was applied to improve the performance. 

The CelebFaces Attributes Dataset (CelebA) is a large-scale face attributes dataset with 202,599 celebrity images, each in color with an original definition of 178 x 218. Each image in the CelebA dataset is in “.jpg” format was imported and converted into 1D vector format. The evaluation will be performed by calculating classification error rate (𝐸 = 𝑁𝑤𝑟𝑜𝑛𝑔/𝑁𝑉) in under the one-hot coding scheme for “Eyeglasses” attribute from list_attr_celeba.txt file indicating whether the person in the picture is wearing glasses or not.

# Project 5

The goal of this study is to classify the operators in a percutaneous nephrolithotomy (PCNL), a minimally-invasive procedure for stone removal from the kidney using a small puncture wound on the skin, simulation into two groups of experts and novices using kinematic features derived from the data of the simulation tool. 

Fourteen participants performed the simulation that included 1 undergraduate student, 4 medical students, 1 post-doctoral research fellow, 5 residents, 2 clinical fellows, and 1 faculty. They were labeled into two levels of minimal/none and residence/training in PCNL experience corresponding to novices and experts, respectively. In a quiet test environment using the two stone kidney model, the player removed both stones to complete the simulation. 10 kinematic features were extracted and normalized from the simulation system data that included the mean value for the left and right tool in cumulative task time, path length, mean and variance velocity, and mean and variance orientation in each of x, y, and z axes. The summary statistics by class for each feature is provided in Table 1 (see the Word document). Linear discriminant analysis with diagonal covariance regularization was employed for the classification. 

The surgeon skill classification was performed on 14 PCNL cases. Linear discriminant model was trained and validated using the aforementioned features by a 5-fold cross-validation with a test accuracy of 85.7% in classifying the skills.

# Project 6-10

The goal of these projects is to use a graphical interface (GUI) to recognize and classify the handwritten digits using different machine learning algorithms, e.g., SVM, Rubine, 1$ Recognizer, Image-based Classifier, and AdaBoost.

# Project 11

The goal of this project is to examine different machine learning methods, e.g., HMM and Rubine, for the classification and recognition of the US and NY zipcodes.

# Project 12 

The attached CSV file lists the customer, date, and dollar value of orders placed at a store in 2017. The gender of each customer is also provided.

The Python code will do the following:

A) Assemble a dataframe with one row per customer and the following columns:
    * customer_id
    * gender
    * most_recent_order_date
    * order_count (number of orders placed by this customer)

B) Plot the count of orders per week.

C) Compute the mean order value for gender 0 and for gender 1 and determine the significance of difference. 


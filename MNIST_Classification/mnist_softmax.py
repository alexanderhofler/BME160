#!/usr/bin/env python3
# Name: Alexander Hoefler (ahoefler)/ Tensorflow
# Group: Alexander Hoefler and Lisa DeFeo (ldefeo) (used tensorflow)


# TO USE: 1. Make sure tensorflow is installed and working on your device
#              -     For installation, see: https://www.tensorflow.org/versions/r1.0/install/
#         2. After installing, copy/paste the code into the command prompt to initialize an interactive
#             tensorflow session and run the program/train the weights.

"""
Uses tensorflow tutorial to create a deep learning convolution layer for image recognition.
Set up to train so that the end result gets ~0.98/1 accuracy for predictions/recignizing features.
Trained using MNIST dataset.
"""

from tensorflow.examples.tutorials.mnist import input_data

mnist = input_data.read_data_sets('MNIST_data', one_hot=True)

import tensorflow as tf

sess = tf.InteractiveSession()

x = tf.placeholder(tf.float32, shape=[None, 784])  # placeholder - size of the pixel matrix (28x28)

y_ = tf.placeholder(tf.float32, shape=[None, 10])  # placeholder - possible output numbers (0-9)

W = tf.Variable(tf.zeros([784, 10]))  # matrix (784,10) variable

b = tf.Variable(tf.zeros([10]))  # vector set for all 10 classes (0-9)

sess.run(tf.global_variables_initializer())

y = tf.matmul(x, W) + b

cross_entropy = tf.reduce_mean(  # takes the average of the softmax sums

    tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))  # applies softmax activation and summates them

train_step = tf.train.GradientDescentOptimizer(0.5).minimize(
    cross_entropy)  # trains the weight with a 0.5 step length - sets gradient

for _ in range(1000):
    batch = mnist.train.next_batch(100)

    train_step.run(feed_dict={x: batch[0], y_: batch[1]})

correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_,
                                                         1))  # returns boolean (y or n) for if the number matches the groundtruth value

accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))  # measures the accuracy of the model

"""
Building the multi-layered convolution network to boost the accuracy.
"""


def weight_variable(shape):
    initial = tf.truncated_normal(shape,
                                  stddev=0.1)  # initializing neurons - non-zero in order to avoid dying ReLu problem
    return tf.Variable(initial)


def bias_variable(shape):
    initial = tf.constant(0.1, shape=shape)
    return tf.Variable(initial)


def conv2d(x, W):
    return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')  # setup step/stride size for convolution


def max_pool_2x2(x):
    return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],  # pools the convoluted filters to reduce size - ignores exact position
                          strides=[1, 2, 2, 1], padding='SAME')


W_conv1 = weight_variable([5, 5, 1, 32])  # set size of the matrix and amount of input/output channels
b_conv1 = bias_variable([32])  # vector for every output channel

x_image = tf.reshape(x,
                     [-1, 28, 28, 1])  # adding new dimensions to x-variable (height, width, color channels -e.g. RGB)

h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)  # activation function (ReLu - half-rectified)
h_pool1 = max_pool_2x2(h_conv1)  # pooling again, 2x2 step/stride

"""
Add a second convolution layer and fully connected layer
"""

W_conv2 = weight_variable([5, 5, 32, 64])
b_conv2 = bias_variable([64])

h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
h_pool2 = max_pool_2x2(h_conv2)

W_fc1 = weight_variable([7 * 7 * 64, 1024])  # reduced and laid out so the model can evaluate the entire picture
b_fc1 = bias_variable([1024])

h_pool2_flat = tf.reshape(h_pool2,
                          [-1, 7 * 7 * 64])  # fully connected layer reshapes into weight matrix multiplied by vectors
h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

keep_prob = tf.placeholder(tf.float32)  # avoid overfitting - neurons are dropped (chosen at "random")
h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)  # "turns off" some weights, to keep other weights from relying on it
#  keeping weights from being too sensitive to certain data sets

W_fc2 = weight_variable([1024, 10])
b_fc2 = bias_variable([10])

y_conv = tf.matmul(h_fc1_drop, W_fc2) + b_fc2

"""
Changing the steep (0.5) gradient descent from above to an ADAM Optimizer - adaptive vs. the old single rate learner (0.5).
Logging accuracy every 100 iterations. Improved/controlled dropout rate. Improved accuracy (99+ %).
"""
cross_entropy = tf.reduce_mean(
    tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y_conv))
train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)
correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
sess.run(tf.global_variables_initializer())
for i in range(20000):
    batch = mnist.train.next_batch(50)
    if i % 100 == 0:
        train_accuracy = accuracy.eval(feed_dict={
            x: batch[0], y_: batch[1], keep_prob: 1.0})
        print("step %d, training accuracy %g" % (i, train_accuracy))
    train_step.run(feed_dict={x: batch[0], y_: batch[1], keep_prob: 0.5})

print("test accuracy %g" % accuracy.eval(feed_dict={
    x: mnist.test.images, y_: mnist.test.labels, keep_prob: 1.0}))

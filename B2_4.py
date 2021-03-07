import scipy.io as sc
import matplotlib.pyplot as plt
import numpy as np

data = sc.loadmat('mnist_sevens_nines.mat')

plt.figure(1)

for i in range(16):
    plt.subplot(4,4,i+1)
    sample = data['X_train'][:,i]
    plt.imshow(sample.reshape(28,28).T, cmap='gray')

#plt.axes('off')
plt.tight_layout()
plt.show()

plt.figure(2)
for i in range(16):
    plt.subplot(4,4,i+1)
    sample = data['X_train'][:,i+900]
    plt.imshow(sample.reshape(28,28).T, cmap='gray')

#plt.axes('off')
plt.tight_layout()
plt.show()


P = data['X_train'].shape[1]
max_itr = 5000
alpha = 0.1/P

n_inp = data['X_train'].shape[0]
n_h1 = 100 #number hidden layer 1 neurons
n_h2 = 80 #^layer 2

weight_scale = 0.01
w3 = weight_scale*np.random.normal(size=(1,n_h2))
w2 = weight_scale*np.random.normal(size=(n_h2,n_h1))
w1 = weight_scale*np.random.normal(size=(n_h1,n_inp))

#NEED TO SHUFFLE TRAINING DATA?
L_train = np.zeros(max_itr)
L_test = np.zeros(max_itr)
acc_train = np.zeros(max_itr)
acc_test = np.zeros(max_itr)

print('Training network...')
for i in range(max_itr):

    z1 = w1 @ data['X_train']
    h1 = np.maximum(z1, 0)
    z2 = w2 @ h1
    h2 = np.maximum(z2, 0)
    yh = w3 @ h2 

    _z1 = w1 @ data['X_test']
    _h1 = np.maximum(_z1, 0)
    _z2 = w2 @ _h1
    _h2 = np.maximum(_z2, 0)
    _yh = w3 @ _h2 

    L_train[i] = 0.5 * np.linalg.norm(data['y_train'] - yh)**2
    acc_train[i] = np.mean(np.sign(yh) == np.sign(data['y_train']))

    L_test[i] = 0.5 * np.linalg.norm(data['y_test'] - _yh)**2
    acc_test[i] = np.mean(np.sign(_yh) == np.sign(data['y_test']))

    if i % 100 == 0:
        print('Epoch:\t', i)
        print('Training loss:\t', L_train[i])
        print('Training acc:\t', acc_train[i])
        print('Test loss:\t', L_test[i])
        print('Test acc:\t', acc_test[i])

    e = data['y_train'] - yh
    d3 = e
    d2 = w3.T @ d3 * (z2 > 0)
    d1 = w2.T @ d2 * (z1 > 0) #check z or h

    g3 = -d3 @ h2.T
    g2 = -d2 @ h1.T
    g1 = -d1 @ data['X_train'].T

    w3 -= alpha*g3
    w2 -= alpha*g2
    w1 -= alpha*g1

print('Done training.')
t = np.arange(max_itr)

plt.figure(3)
plt.subplot(121)
plt.plot(t, L_train/P, L_test/data['X_test'].shape[1], linewidth=2)
plt.legend(['Train loss', 'Test loss'])
plt.ylabel('MSE')
plt.subplot(122)
plt.plot(t, acc_train, acc_test, color='r', linewidth=2)
plt.legend(['Train accuracy', 'Test accuracy'])
plt.xlabel('Epochs')
plt.ylabel('Accuracy %')
plt.show()




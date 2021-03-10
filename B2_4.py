import scipy.io as sc
import numpy as np
import csv
import matplotlib.pyplot as plt

def runNet(data, n_h1=100, n_h2=80, weight_scale=0.01, save=True):
    n_inp = data['X_train'].shape[0]

    w3 = weight_scale*np.random.normal(size=(1,n_h2))
    w2 = weight_scale*np.random.normal(size=(n_h2,n_h1))
    w1 = weight_scale*np.random.normal(size=(n_h1,n_inp))

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
    t = np.arange(1,max_itr+1)

    if save:
        print('Saving data...')
        data_out = zip(t, L_train/P, L_test/data['X_test'].shape[1], acc_train, acc_test)
        filename = 'logs/B2_4_data_w{x}.csv'.format(x = weight_scale)
        np.savetxt(filename, [p for p in data_out], delimiter=',', fmt='%f')
        print('Saved.')

    return data_out

#LOAD DATA

data = sc.loadmat('data/mnist_sevens_nines.mat')

#DATA VISUALISATION

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


#WEIGHTSCALE 0.01 EXPERIMENT

runNet(data)

##WEIGHTSCALE 0.1 EXPERIMENT

runNet(data, weight_scale=0.1)




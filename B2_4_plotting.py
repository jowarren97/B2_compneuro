import matplotlib.pyplot as plt
import numpy as np

data1 = np.genfromtxt('logs/B2_4_data_w0.01.csv',delimiter=',')
data2 = np.genfromtxt('logs/B2_4_data_w0.1.csv', delimiter=',')

def plotPerf(data, weightscale, save=True, plot=True):
    fig, ax = plt.subplots(1,2, figsize=(10,6))
    #plt.figure(figsize=(6,6))
    #plt.subplot(211)
    ax[0].plot(data[:,0], data[:,1], color='r', linewidth=1) #t, L_train/P
    ax[0].plot(data[:,0], data[:,2], color='b', linewidth=1) #t, L_test/X_train.shape
    ax[0].legend(['Train loss', 'Test loss'])
    ax[0].set_xlabel('Epochs trained')
    ax[0].set_ylabel('Loss (MSE)')
    #ax[0].get_xaxis().set_visible(False)
    ax[0].set_xlim([0,data[-1,0]])
    ax[0].set_ylim(bottom=0)

    ax[1].plot(data[:,0], data[:,3], color='r', linewidth=1) #t, acc_train
    ax[1].plot(data[:,0], data[:,4], color='b', linewidth=1) #t, acc_test
    ax[1].legend(['Train accuracy', 'Test accuracy'])
    ax[1].set_xlabel('Epochs trained')
    ax[1].set_ylabel('Accuracy')
    ax[1].set_xlim([0,data[-1,0]])
    ax[1].set_ylim(bottom=0.5) #chance

    fig.suptitle("Weightscale = {x}".format(x=weightscale))

    if save:
        plt.savefig("figs/png/B2_4_w{x}.png".format(x=weightscale))
    if plot:
        plt.show()


plotPerf(data1, 0.01)
plotPerf(data2, 0.1)

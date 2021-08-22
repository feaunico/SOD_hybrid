
import numpy as np


import matplotlib
import matplotlib.pyplot as plt


def plotPCA(fichier):
    print fichier

    fx = open(fichier)
    ct = fx.readlines()
    fx.close()


    color = {'NA1':'b','NA2':'g','EU1':'k','EU2':'orange','unknown_Guillaume':'r'}
    perc = ct[1].split('\t')[6].split(' ')




    fig, ax = plt.subplots(figsize=(9, 6))

    ax0 = plt.subplot2grid((2, 2), (0, 0), rowspan=1)
    ax1 = plt.subplot2grid((2, 2), (0, 1), rowspan=1)
    ax2 = plt.subplot2grid((2, 2), (1, 0), rowspan=1)
    ax3 = plt.subplot2grid((2, 2), (1, 1), rowspan=1)

    for x in ct:
        ax0.plot(float(x.split('\t')[1]),float(x.split('\t')[2]),marker = '.', color = color[x.split('\t')[4].replace('\n','4')])
        ax1.plot(float(x.split('\t')[1]), float(x.split('\t')[3]), marker='.',
                 color=color[x.split('\t')[4].replace('\n', '4')])
        ax2.plot(float(x.split('\t')[2]), float(x.split('\t')[3]), marker='.',
                 color=color[x.split('\t')[4].replace('\n', '4')])

    ax3.set_xlim(0,10)
    ax3.set_ylim(0,10)
    ax3.plot(1,7,marker = '.', color = color['NA1'])
    ax3.text(1.5,7,'NA1',va = 'center')

    ax3.plot(1,5.5,marker = '.', color = color['NA2'])
    ax3.text(1.5,5.5,'NA2',va = 'center')

    ax3.plot(4,7,marker = '.', color = color['EU1'])
    ax3.text(4.5,7,'EU1',va = 'center')

    ax3.plot(4,5.5,marker = '.', color = color['EU2'])
    ax3.text(4.5,5.5,'EU2',va = 'center')

    ax3.plot(1,4,marker = '.', color = color['unknown_Guillaume'])
    ax3.text(1.5,4,'Phra1_7964 and Phra1_7974',va = 'center')
    ax3.axis('off')
    ax0.set_yticks([0])
    ax0.set_xticks([0])
    ax1.set_yticks([0])
    ax1.set_xticks([0])
    ax2.set_yticks([0])
    ax2.set_xticks([0])

    ax0.set_xticklabels(['PC1 (' + perc[0] + '%)'],fontsize=8)
    ax0.set_yticklabels(['PC2\n(' + perc[1] + '%)'],fontsize=8)

    ax1.set_xticklabels(['PC1 (' + perc[0] + '%)'],fontsize=8)
    ax1.set_yticklabels(['PC3\n(' + perc[2] + '%)'],fontsize=8)
    ax2.set_xticklabels(['PC2 (' + perc[1] + '%)'],fontsize=8)
    ax2.set_yticklabels(['PC3\n(' + perc[2] + '%)'],fontsize=8)


    ax0.set_ylim(ax0.get_ylim())
    ax0.set_xlim(ax0.get_xlim())
    ax1.set_ylim(ax1.get_ylim())
    ax1.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax2.get_ylim())
    ax2.set_xlim(ax2.get_xlim())

    ax0.plot((0,0),ax0.get_ylim(),'k--',linewidth=0.7)
    ax0.plot(ax0.get_xlim(),(0, 0), 'k--',linewidth=0.7)
    ax1.plot((0,0),ax1.get_ylim(),'k--',linewidth=0.7)
    ax1.plot(ax1.get_xlim(),(0, 0), 'k--',linewidth=0.7)
    ax2.plot((0,0),ax2.get_ylim(),'k--',linewidth=0.7)
    ax2.plot(ax2.get_xlim(),(0, 0), 'k--',linewidth=0.7)


    plt.savefig(fichier.replace('txt','png'),dpi=800,  format='png')


plotPCA('NoPrunning.txt')
plotPCA('LD02.txt')
plotPCA('LD05.txt')
plotPCA('LD07.txt')
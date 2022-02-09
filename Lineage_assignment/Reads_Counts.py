import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

fx = open('PR_hyb.recode.vcf')
ct = fx.readlines()
fx.close()

chi1, chi2 = [], []
n, t = 0, 0
nn, tt = 0, 0
ls1, ls2 = [], []
for x in ct:
    if x[0] != '#':
        genot1, genot2 = x.split('\t')[71], x.split('\t')[72]
        if genot1.split(':')[0] == '0/1' or genot1.split(':')[0] == '0|1' or genot1.split(':')[0] == '1|0':
            if int(genot1.split(':')[3]) >=99:
                t = t + 1
                val1, val2 = int(genot1.split(':')[1].split(',')[0]), int(genot1.split(':')[1].split(',')[1])
                if val1 <= val2:
                    vv = float(val1) / (val1 + val2)
                    ls1.append(vv)
                else:
                    vv = float(val2) / (val1 + val2)
                    ls1.append(vv)
                exp2 = 0.5 * (val1 + val2)
                exp1 = exp2
                p = stats.chisquare([val1,val2], [exp1,exp2])
                if p[1] > 0.05:
                    n = n + 1
                chi1.append(p[1])

        if genot2.split(':')[0] == '0/1' or genot2.split(':')[0] == '0|1' or genot2.split(':')[0] == '1|0':
            if int(genot1.split(':')[3]) >=99:
                tt = tt + 1
                val1, val2 = int(genot2.split(':')[1].split(',')[0]), int(genot2.split(':')[1].split(',')[1])
                if val1 <= val2:
                    vv = float(val1) / (val1 + val2)
                    ls2.append(vv)
                else:
                    vv = float(val2) / (val1 + val2)
                    ls2.append(vv)
                exp2 = 0.5 * (val1 + val2)
                exp1 = exp2
                p = stats.chisquare([val1,val2], [exp1,exp2])
                if p[1] > 0.05:
                    nn = nn + 1
                chi2.append(p[1])



print np.mean(ls1), np.std(ls1)
print np.mean(ls2), np.std(ls2)
print n, t, (float(n)/t) *100
perc1 = round((float(n)/t) *100,2)
print nn, tt, (float(nn)/tt) *100
perc2 = round((float(nn)/tt) *100,2)


############################### BUILD FIGURE S1 ##################################

fig, ax = plt.subplots(figsize=(7, 4))
ax0 = plt.subplot2grid((1, 2), (0, 0), colspan=1, rowspan=1)
ax1 = plt.subplot2grid((1, 2), (0, 1), colspan=1, rowspan=1)

ax0.hist(ls1, color = 'grey', bins = 50, )
ax1.hist(ls2, color = 'grey', bins = 50, )


ax0.set_xticks([0,0.25, 0.5])
ax1.set_xticks([0,0.25,0.5])
ax0.set_title('16-237-021',fontsize = 9)
ax1.set_title('16-284-032',fontsize = 9)
ax0.set_ylabel('Number of SNP loci')
ax0.set_xlabel('Proportion of reads\n for minor allele')
ax1.set_xlabel('Proportion of reads\n for minor allele')


ax0.text(0.025,11000,str(perc1) + '% heterozygous loci\nwith Chi-square prob > 0.05',fontsize = 8)
ax1.text(0.025,15000,str(perc2) + '% heterozygous loci\nwith Chi-square prob > 0.05',fontsize = 8)

plt.tight_layout()
plt.savefig('Figure_S1.png',dpi=800,  format='png')
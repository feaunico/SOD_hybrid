
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sklearn.naive_bayes import MultinomialNB
from sklearn.naive_bayes import GaussianNB
import os, sys
import time



def cluster(data, maxgap):
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups




start_time = time.time()
print len(sys.argv), sys.argv

if len(sys.argv) == 5:
    reads = sys.argv[3] + ' ' + sys.argv[4]
    training_set = sys.argv[2]
    ref = sys.argv[1]
elif len(sys.argv) == 4:
    reads = sys.argv[3]
    training_set = sys.argv[2]
    ref = sys.argv[1]
else:
    print 'usage: python2.7 ASSIGM.py <ref_genome> <training_set> <fastq_reads_1> <fastq_reads_2>'
    print
    print 'Convert SRA data in fastq using => fastq-dump --split-files name_of_SRA_file '
    exit()

print ref
print training_set
print reads



fx = open(ref)
ct = fx.read().split('>')
fx.close()
ct = ct[1:]

dico = {}
for x in ct:

    dico[x.split('\n',1)[0]] = x.split('\n',1)[1].replace('\n','').replace('\r','') + 'N'*500




########################################################
##                                                    ##
#       Training of the naive bayesian classifier      #
#                with the input SNP dataset            #
##                                                    ##
########################################################

fx = open(training_set)
snp = fx.readlines()
fx.close()

header = snp[0].split('\t')[2:-1] + [snp[0].split('\t')[-1].replace('\n','').replace('\r','')]
loci = []

snp = snp[1:]
snp.sort()

matrix = []
for x in snp:
    if x.split('\t')[0] != '\n' and x.split('\t')[0] != '':
        ls = []
        Line = x.split('\t')
        Line = Line[2:]

        for y in Line:
            y = y.replace('\n','').replace('\r','')
            y = y.replace('0|0','1').replace('0|1','2').replace('1|1','3').replace('2|2','3').replace('3|3','3').replace('4|4','3').replace('5|5','3').replace('6|6','3').replace('7|7','3').replace('8|8','3').replace('9|9','3').replace('.|.','1')
            ls.append(int(y))
        matrix.append(ls)

        loci.append(x.split('\t')[0] + '|' + x.split('\t')[1])

print len(matrix), len(matrix[0])

matrix = np.array(matrix)
matrix = matrix.transpose()

clf = MultinomialNB()
#clf = GaussianNB()
clf.fit(matrix, header)
print header

uniq = []
for x in header:
    if x not in uniq:
        uniq.append(x)


########################################################
##                                                    ##
#     Generation of the synthetic reference genome     #
##                                                    ##
########################################################

scaffolds = []
for x in loci:
    if x.split('|')[0] not in scaffolds and x.split('|')[0] != '\n':
        scaffolds.append(x.split('|')[0])
scaffolds.sort()


coord = []
seq = ''
i = 500
corresp = {}


for sca in scaffolds:
    ls = []
    for y in loci:
        if y.split('|')[0] == sca:
            ls.append(int(y.split('|')[1]))
    ls.sort()
    groups = cluster(ls, maxgap=500)

    for g in groups:

        if len(g) == 1:

            seq = seq + dico[sca][g[0] - 500:g[0] + 500]
            coord.append(i-1)
            corresp[str(i-1)] = sca + '|' + str(g[0])
        else:
            seq = seq + dico[sca][g[0] - 500:g[-1] + 500]
            coord.append(i-1)
            corresp[str(i - 1)] = sca + '|' + str(g[0])


            j = 1
            while j < len(g):
                nb = g[j] - g[j-1]
                k = i + nb -1
                coord.append(k)
                corresp[str(k)] = sca + '|' + str(g[j])
                i = i + nb
                j = j + 1
        seq = seq + 'N' * 1000
        i = i + 2000

fout = open('Synthetic_genome.fas','w')
fout.write('>Synth\n' + seq + '\n')
fout.close()


interm_time = time.time() - start_time

####################################################################
##                                                                ##
#     Alignment of the reads on the synthetic reference genome     #
##                                                                ##
####################################################################


os.system('samtools faidx Synthetic_genome.fas')
os.system('bwa index Synthetic_genome.fas')
os.system('bwa mem -T 24 Synthetic_genome.fas ' + reads + ' > Rs.sam')

os.system('samtools view -b -F 4 Rs.sam > Rs_mapped.bam')
os.system('samtools sort Rs_mapped.bam -o Rs_mapped_sorted.bam')
os.system('samtools index Rs_mapped_sorted.bam')

outls = []
logls = []

wr = 0
averageDP, averageMQ = [], []
for x in corresp:
    x = int(x)
    print x
    res = os.popen('bcftools mpileup -Ou -r Synth:' + str(x+1) + '-' + str(x+1) + ' -f Synthetic_genome.fas Rs_mapped_sorted.bam | bcftools call -m').readlines()
    print res[-1]
    try:
        MQ = res[-1].split('MQ=')[1].split('\t')[0]
        if MQ == '.':
            MQ = 0.0
        dicoVal = {'DP':float(res[-1].split('DP=')[1].split(';')[0]), 'MQ': MQ, 'MQ0F': float(res[-1].split('MQ0F=')[1].split(';')[0])}
        genot = res[-1].split('\t')[-1].replace('\n','')
        genot = genot.split(':')[0]


        averageDP.append(dicoVal['DP'])
        averageMQ.append(float(MQ))

    except:
        genot = './.'
        pass
    print 'GENOT', genot
    if genot == 'Rs_mapped_sorted.bam' or genot == './.':
        genot = '0/0'
        wr = wr + 1
    elif genot =='1/2':
        genot = '0/1'
    print '**', genot


    outls.append(corresp[str(x)] + '|' + genot)
    try:
        logls.append((int(x), str(corresp[str(x)]) + '\t' + str(genot) + '\t' + str(dicoVal['DP']) + '\t' + str(dicoVal['MQ']) + '\t' + str(dicoVal['MQ0F']) + '\n'))
    except:
        logls.append((int(x), str(corresp[str(x)]) + '\t' + str(genot) + '\t0\t0\tNAN\n'))

logls.sort(key=lambda a: a[0])

flog = open('Log.txt','w')

flog.write('Location\tRecalibrated_Loc\tGenotype\tDepth\tMapping Quality\tFraction of MQ0 reads\n')
for x in logls:
    flog.write(x[1])

flog.close()
outls.sort()
print outls


ls = []
for x in outls:
    genot = x.split('|')[-1].replace('0/0','1').replace('0/1','2').replace('1/1','3').replace('\n','').replace('\r','')

    print 'voici genot', x, genot
    ls.append(int(genot))

pred = clf.predict([ls])
pred1 = clf.predict_proba([ls])
print pred
print pred1

fout = open('log.assign','w')
fout.write(str(pred) + '\t')
fout.write(str(clf.classes_) + '\t')
fout.write(str(pred1) + '\n')
fout.close()

print 'Average depth', np.mean(averageDP), '+/-', np.std(averageDP)
print 'Average MQ', np.mean(averageMQ), '+/-', np.std(averageMQ)


print len(coord) - wr, 'SNPs used out of', len(coord)

############################### BUILD FIGURE S1 ##################################

fig, ax = plt.subplots(figsize=(7, 4))
ax0 = plt.subplot2grid((1, 2), (0, 0), colspan=1, rowspan=1)
ax1 = plt.subplot2grid((1, 2), (0, 1), colspan=1, rowspan=1)
ax0.hist(averageDP, color = 'grey', bins = 50, )
ax1.hist(averageMQ, color = 'grey', bins = 50, )
ax0.set_ylabel('Number of SNP loci')
ax0.set_xlabel('Depth')
ax1.set_xlabel('Mapping quality')
plt.tight_layout()
plt.savefig('Figure_S1.png',dpi=800,  format='pdf')
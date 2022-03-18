


import os
import numpy as np


def clean(fichier,seqname):
    fto = open(seqname,'w')
    ls = []
    naMe = []
    for x in fichier:
        ls.append(list(x.split('\n',1)[1].replace('\n','')))
        naMe.append(x.split('\n',1)[0])


    lss = [list(i) for i in zip(*ls)]

    nls =[]
    for x in lss:

        locus = []
        uniq = []

        for y in x:
            if y not in uniq:
                uniq.append(y)

        if len(uniq) >= 2:

            if 'N' not in uniq:
                print uniq
                nls.append(x)


    print 'nls', nls[0]
    lss = [list(i) for i in zip(*nls)]
    print '==>', lss[0]
    eu1,eu2,na1,na2 = [], [], [], []
    for x,y in zip(lss,naMe):
        print y
        if y.split('_')[-2] == 'EU1':
            if "".join(x) not in eu1:
                eu1.append("".join(x))
                fto.write(y + '\n' + "".join(x)+ '\n')
        elif y.split('_')[-2] == 'NA1':
            if "".join(x) not in na1:
                na1.append("".join(x))
                fto.write(y + '\n' + "".join(x)+ '\n')
        elif y.split('_')[-2] == 'EU2':
            if "".join(x) not in eu2:
                eu2.append("".join(x))
                fto.write(y + '\n' + "".join(x)+ '\n')
        elif y.split('_')[-2] == 'NA2':
            if "".join(x) not in na2:
                na2.append("".join(x))
                fto.write(y + '\n' + "".join(x)+ '\n')
        else:
            fto.write(y + '\n' + "".join(x) + '\n')
    fto.close()


def parse(liste ,coord, head , gene):


    n = 0
    for y in liste:
        if '|' in y.split('\t')[71].split(':')[0]:
            n = n + 1


    if n >= 15 and n <20:
        NA1, NA2, EU1, EU2 = [], [], [], []
        fout = open(gene + '.fas', 'w')
        lfl = []
        print n

        seq = fas[coord[0]]

        i = 0
        head = head[9:]
        while i < len(head):
            name = '>' + head[i].replace('\n','') + '_' + names[head[i].replace('\n','')]

            hap1, hap2 = list(seq), list(seq)

            for ent in liste:

                pos = int(ent.split('\t')[1])
                alt, ref = ent.split('\t')[4], ent.split('\t')[3]
                ent = ent.split('\t')[9:]

                genot = ent[i]
                if '.' in genot.split(':')[0]:
                    hap1[pos] = 'N'
                    hap2[pos] = 'N'

                elif genot.split(':')[0] == '0/1':
                    hap1[pos] = 'N' #deg[alt+ref]
                    hap2[pos] = 'N' #deg[alt+ref]
                elif genot.split(':')[0] == '1|1' or genot.split(':')[0] == '1/1':
                    hap1[pos] = alt
                    hap2[pos] = alt
                elif genot.split(':')[0] == '0|0' or genot.split(':')[0] == '0/0':
                    hap1[pos] = ref
                    hap2[pos] = ref
                elif genot.split(':')[0] == '0|1':
                    hap1[pos] = ref
                    hap2[pos] = alt

                elif genot.split(':')[0] == '1|0':
                    hap1[pos] = alt
                    hap2[pos] = ref

            hap_1 = name + '_1\n' + "".join(hap1)[coord[1]:coord[2]] + '\n'
            hap_2 = name + '_2\n' + "".join(hap2)[coord[1]:coord[2]] + '\n'



            if names[head[i].replace('\n','')] == 'NA1':

                if "".join(hap1)[coord[1]:coord[2]] not in NA1:
                    NA1.append("".join(hap1)[coord[1]:coord[2]])
                    fout.write(hap_1)
                    lfl.append(hap_1)
                if "".join(hap2)[coord[1]:coord[2]] not in NA1:
                    NA1.append("".join(hap2)[coord[1]:coord[2]])
                    fout.write(hap_2)
                    lfl.append(hap_2)
            elif names[head[i].replace('\n','')] == 'NA2':
                if "".join(hap1)[coord[1]:coord[2]] not in NA2:
                    NA2.append("".join(hap1)[coord[1]:coord[2]])
                    fout.write(hap_1)
                    lfl.append(hap_1)
                    lfl.append(hap_2)
                if "".join(hap2)[coord[1]:coord[2]] not in NA2:
                    NA2.append("".join(hap2)[coord[1]:coord[2]])
                    fout.write(hap_2)
                    lfl.append(hap_2)
            elif names[head[i].replace('\n','')] == 'EU1':
                if "".join(hap1)[coord[1]:coord[2]] not in EU1:
                    EU1.append("".join(hap1)[coord[1]:coord[2]])
                    fout.write(hap_1)
                    lfl.append(hap_1)
                if "".join(hap2)[coord[1]:coord[2]] not in EU1:
                    EU1.append("".join(hap2)[coord[1]:coord[2]])
                    fout.write(hap_2)

                    lfl.append(hap_2)



            elif names[head[i].replace('\n','')] == 'EU2':
                if "".join(hap1)[coord[1]:coord[2]] not in EU2:
                    EU2.append("".join(hap1)[coord[1]:coord[2]])
                    fout.write(hap_1)
                    lfl.append(hap_1)


                if "".join(hap2)[coord[1]:coord[2]] not in EU2:
                    EU2.append("".join(hap2)[coord[1]:coord[2]])
                    fout.write(hap_2)

                    lfl.append(hap_2)


            else:
                fout.write(hap_1)
                fout.write(hap_2)
                lfl.append(hap_1)
                lfl.append(hap_2)


            i = i + 1
        fout.close()

        clean(lfl, gene + '_clean.fas')

        os.system('raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s ' + gene + '.fas -n ' + gene + '.tree')
        os.system('raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s ' + gene + '_clean.fas -n ' + gene + '_clean.tree')




deg = {'AC':'M','CA':'M','AG':'R','GA':'R','CT':'Y','TC':'Y','GC':'S','CG':'S','AT':'W','TA':'W','GT':'K','TG':'K'}



fx = open('PHRA102.gff')
ct = fx.readlines()
fx.close()

gff = []
for x in ct[1:]:
    if x.split('\t')[2] == 'gene':
        gff.append((x.split('\t')[0],int(x.split('\t')[3]),int(x.split('\t')[4]), x.split('\t')[8].replace('ID=','').replace('\n','') ))
print 'gff done'



fx = open('PR-102_v3.1.fasta')
ct = fx.read().split('>')
fx.close()

fas = {}
for x in ct[1:]:
    fas[x.split('\n')[0]] = x.split('\n',1)[1].replace('\n','').replace('\r','')
print 'fasta done'


fx = open('listeofisolates')
ct = fx.readlines()
fx.close()

names = {}
for x in ct:
    names[x.split('\t')[0]] = x.split('\t')[1].replace('\n','')
print 'isolate names done'


fx = open('PR_hyb.recode.vcf')
vcf = fx.readlines()
fx.close()

for x in vcf[78:82]:
    if x.split('\t')[0] == '#CHROM':
        tete = x.split('\t')

print 'starting gff'

m = 0
vcf = vcf[80:]
for x in gff:
    #print x
    ls = []

    i = m

    while i < len(vcf):
        line = vcf[i]
        m = m + 1
        if line.split('\t')[0] == x[0] and int(line.split('\t')[1]) >= x[1] and int(line.split('\t')[1]) <= x[2]:
                ls.append(line)


        if int(line.split('\t')[1]) > x[2]:
            if len(ls) >= 5:
                parse(ls,x,tete,x[-1])
            i = len(vcf)
        i = i + 1







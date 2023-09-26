from os import system
versionfile=open('BlasLike/Licence.cs')
newfile=open('Licence.cs','w')
while 1:
    line=versionfile.readline()
    if len(line)==0:break
    tt=line.find('int revision')
    if line.find('int revision')!= -1:
        revision=int(line.split('=')[1].replace(';',''))
        print('Old revision %d'%revision)
        revision+=1
        print('New revision %d'%revision)
        line=('\t\tpublic readonly int revision=%d;\n' % (revision))
    newfile.write(line)
newfile.close()
versionfile.close()
system('\\cygwin\\bin\\mv Licence.cs BlasLike\\Licence.cs');


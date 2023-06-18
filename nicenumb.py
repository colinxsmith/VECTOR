from math import *

def betterdelt(delt):
    if delt>=5:delt=5
    elif delt>=2:delt=2
    elif delt>=1:delt=1
    return delt

def nicenumb(a):
    delt=abs(a)
    mul=1
    if delt>1:
        while delt>10:
            delt/=10
            mul-=1
        delt=round(delt)
        delt=betterdelt(delt)
        while mul<1:
            delt*=10
            mul+=1
        return delt
    else:
        while delt<1:
            delt*=10
            mul+=1
        delt=round(delt)
        delt=betterdelt(delt)
        while mul>1:
            delt/=10
            mul-=1
        return delt
k=-0.0857326569979281
xmin=k
xmax=0
print(k,nicenumb(xmax-xmin))
dx=nicenumb(xmax-xmin)
tick=xmin/dx
print(tick,round(tick-.5))
xmin=round(tick-0.5)*dx
print(nicenumb(xmax-xmin),xmin,k,xmax)
k=34.34567678
xmin=0
xmax=k
print(k,nicenumb(xmax-xmin))
dx=nicenumb(xmax-xmin)
tick=xmax/dx
print(tick,round(tick+0.5))
xmax=round(tick+0.5)*dx
print(nicenumb(xmax-xmin),xmin,k,xmax)
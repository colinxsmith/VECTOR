from math import *
def nicenumb(a,ticks=4):
    delt=abs(a/ticks)
    mul=1
    if delt>1:
        while delt>ticks:
            delt/=ticks
            mul-=1
        delt=round(delt)
        while mul<1:
            delt*=10
            mul+=1
        return delt
    else:
        while delt<1:
            delt*=ticks
            mul+=1
        delt=round(delt)
        while mul>1:
            delt/=10
            mul-=1
        return delt
k=-0.0857326569979281
print(k,nicenumb(k))
for d in [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]:
    print('%f %f'% (d*nicenumb(k),k))

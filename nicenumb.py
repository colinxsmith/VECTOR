from math import *
def nicenumb(a):
    delt=abs(a/10)
    mul=1
    if delt>1:
        while delt>10:
            delt/=10
            mul-=1
        delt=round(delt)
        while mul<1:
            delt*=10
            mul+=1
        return delt
    else:
        while delt<1:
            delt*=10
            mul+=1
        delt=round(delt)
        while mul>1:
            delt/=10
            mul-=1
        return delt
for i in [1.23456789e3,1.23456789e2,1.23456789e1,1.23456789,1.23456789e-1,1.23456789e-2,1.23456789e-3,1.23456789e-4,1.23456789e-5]:
    print ('%10f\t%f'%(i,nicenumb(i)))
for i in [-1.23456789e3,-1.23456789e2,-1.23456789e1,-1.23456789,-1.23456789e-1,-1.23456789e-2,-1.23456789e-3,-1.23456789e-4,-1.23456789e-5]:
    print ('%f\t%f'%(i,nicenumb(i)))

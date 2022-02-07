from Optn import *
from re import search
from sys import argv


def testmul(nn, n1, n2, n3, HH, x, y):
    # print(nn)
    H = Hhere
    n = ntrue
    ij = 0
    if len(H) > 0:
        Sym_mult(n, H, x, y)
    for i in range(nn-n):
        y[n+i] = 0


def utility(c, Q, x, tt):
    n = len(x)
    imp = [0]*ntrue
    tt(ntrue, 1, 1, 1, Q, x, imp)
    return dot(c, x)+0.5*dot(imp, x)


class FDATA:
    def __init__(self, ff):
        self.QQ = []
        self.WW = []
        self.bench = []
        keyw = 0
        lastkey = 0
        for line in ff.readlines():
            if line.find('--------------') == 0:
                print('\033[1;1;32mBreak at %s \033[0;m' % line)
                break
            if not keyw:
                keyw = line.strip()
                # This allows numerical data across more than 1 line
                if search('^[0-9-]', keyw[0]):
                    print('\033[1;1;32mExtra line beginning with number for ' +
                          '\033[1;1;31m' + lastkey + '\033[1;1;36m ' + keyw.split()[0] + '\033[0;m')
                    up = getattr(self, lastkey)
                    addon = keyw.split()
                    try:
                        addon = [int(i) for i in addon]
                    except:
                        try:
                            addon = [float(i) for i in addon]
                        except:
                            pass
                    try:up += addon
                    except:
                        up=[up]
                        up+=addon
                    setattr(self, lastkey, up)
                    keyw = 0
            else:
                dd = line.strip().split()
                try:
                    dd = [int(i) for i in dd]
                except:
                    try:
                        dd = [float(i) for i in dd]
                    except:
                        pass
                if len(dd) == 1:
                    setattr(self, keyw, dd[0])
                else:
                    setattr(self, keyw, dd)
                lastkey = keyw
                keyw = 0


if len(argv) > 1:
    Opt = FDATA(open(argv[1]))
else:
    Opt = FDATA(open('./costtestdata'))

A = Opt.AA
C = Opt.CC
L = Opt.LL
U = Opt.UU
n = Opt.n
ntrue = Opt.ntrue
if len(Opt.QQ) > 0 and len(Opt.bench) > 0:
    cextra = []
    Sym_mult(ntrue, Opt.QQ, Opt.bench, cextra)
    for i in range(len(cextra)):
        Opt.CC[i] -= cextra[i]

m = Opt.m
pp = [1, 2, 3]  # Opt.QQ#+[0.0]*int(n*(n+1)-ntrue*(ntrue+1)/2)
print(pp[-1])
print('pp', n*(n+1)/2, len(pp))
print('Q', ntrue*(ntrue+1)/2, len(Opt.QQ))
print('C', n, len(C))
print('bench', ntrue, len(Opt.bench))
print('L', n+m, len(U))
print('U', n+m, len(U))
print('A', n*m, len(A))
LAMBDA = [0]*(n+m)
w = Opt.WW
print('w', n, len(w))
Hhere = Opt.QQ+[0.0]*int(n*(n+1)-ntrue*(ntrue+1)/2)
back = OptAdvanced(n, m, w, A, L, U, C, LAMBDA, testmul, Hhere)
print(back)
print(utility(C, Opt.QQ, w, testmul))

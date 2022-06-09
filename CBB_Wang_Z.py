import numpy as np
import time


#set Adjacency matrix
def Adj(E,V):
    adj = np.zeros((len(V),len(V))) 
    for (i,j) in E:
        adj[i-1][j-1] = 1 
        adj[j-1][i-1] = 1 
    np.fill_diagonal(adj, 0)
    return adj


# $\Gamma(v)$



# Find the v's adjacency vertices set
def gamma(E,V,v):
    adj = Adj(E,V)
    GAMMA = []
    for j in V:
        if adj[v-1][j-1] == 1:
            GAMMA.append(j)
    return GAMMA


# $deg(v)$



# Find the v's adjacency vertices set cardinality
def deg(E,V,v):
    GAMMA = gamma(E,V,v)
    return (len(GAMMA))


# $deg_{V}(v)$



# Find the v's adjacency vertices set cardinality in set V'
def degv(E,V,v,p):
    adj = Adj(E,V)
    r = list(set(V).difference(set(v)))
    for i in r:
        for j in V:
            adj[i-1][j-1] = 0
            adj[j-1][i-1] = 0
    GAMMA = []
    for j in V:
        if adj[p-1][j-1] == 1:
            GAMMA.append(j)
    return len(GAMMA)


# $\Delta(G)$



# find the G's max degree
def delta(E,V):
    adj = Adj(E,V)
    list = np.zeros(len(V))
    for i in V:
        for j in V:
            if adj[i-1][j-1] == 1:
                list[i-1] = list[i-1]+1
    
    return (max(list))


# R_{min} = set of vertices with the minimum degree in R


def RMIN(adj,E,V,R):
    Rmin = np.zeros(len(V))
    mini = []
    for i in R:
        for j in V:
            if adj[i-1][j-1] == 1:
                Rmin[i-1] = Rmin[i-1]+1
    np.place(Rmin,Rmin == 0,float('inf'))
    for i in R:
        if Rmin[i-1] == min(Rmin):
            mini.append(i)         
           
    return mini


# $ex-deg(q)$


def ex_deg(E,V,q):
    exdeg = 0
    for r in gamma(E,V,q):
        exdeg += deg(E,V,r)
    return exdeg


# $Min\{ex-deg(q)|q\in R_{min}\}$




def min_ex_deg(E,V,Rmin):
    list = []
    for q in Rmin:
        list.append(ex_deg(E,V,q))
    Min = min(list)
    for p in Rmin:
        if ex_deg(E,V,p) == Min:
            return p  





def Expand(E,V,R,No):

    NO = No
    R = list(V)
    global Q,Q_max
    while len(R) != 0:   
        p = R[-1]    
        if len(Q)+NO[p] > len(Q_max):
            Q.append(p)
            Rp = [val for val in R if val in gamma(E,V,p)]
            if len(Rp) != 0:
                NUMBER_SORT(E,V,Rp,NO)
                Expand(E,V,Rp,NO)
            elif len(Q) > len(Q_max):
                Q_max = Q
            Q.remove(p)
        else:
            return
        R.remove(p)
        





def NUMBER_SORT(E,V,R,No):
    maxno = 0
    C = {i:[] for i in V}
    while len(R) != 0:
        p = R[0]
        k = 1
        while len([val for val in C[k] if val in gamma(E,V,p)]) != 0:
            k = k+1
        if k > maxno:
            maxno = k
            C[maxno] = [] 
        No[p] = k
        C[k] = list(set([p] + C[k]))
        R.remove(p)
    i = 1
    for k in range(1,maxno+1):
        for j in range(0,len(C[k])):
            R.append(C[k][j])
            i = i+1
    return R,No





def SORT (E,V):
    i = len(V)-1
    adj = Adj(E,V)
    R = list(V)
    VV = {}
    Rmin = RMIN(adj,E,V,R)
    Deg = []
    for j in V:
        Deg.append(deg(E,V,j))
    while len(Rmin) != len(R):
        if len(Rmin) >= 2:
            p = min_ex_deg(E,V,Rmin)
        else:
            p = Rmin[0]
        VV[i] = p
        R.remove(p)
        i = i - 1
        for j in range(1,(len(R))):
            if adj[R[j]-1][p-1] == 1:
                Deg[R[j]-1] = Deg[R[j]-1]-1
        Rmin = RMIN(adj,E,V,R)
    return Rmin, VV





def Start(E,V,VV,No,Rmin):   
    for q in Rmin:
        if degv(E,V,Rmin,q) == len(Rmin)-1:
            Q_max = Rmin
    Expand(E,V,VV,No)
    return Q_max





def NUMBER (E,V, Rmin, VV):

    m = max([No[q] for q in Rmin])
    mmax = len(Rmin) + delta(E,V) - m

    m = m+1
    i = len(Rmin)
    while i <= mmax-1:
        if i > len(VV)-1:
            Q_max = Start(E,V,VV,No,Rmin)

        No[VV[i]] = m
        m = m+1
        i = i+1
    for i in range(int(mmax + 1),len(VV)):
        No[VV[i]] = delta(E,V) + 1
    return Q_max


# C125.9.txt




No = {}
Q = []
Q_max = []
E = np.loadtxt("C125.9.txt" ,dtype=int) # read txt, convert to edge set
V = set(E.flatten()) # get the vertices from edge
start = time.time()
Rmin,VV = SORT(E,V)
Rmin,No = NUMBER_SORT(E,V,Rmin,No)
for i in range(1,len(Rmin)):
    VV[i] = Rmin[i]
Q_max = NUMBER (E,V, Rmin, VV)
end = time.time()
print('the max-clique cardinality is:',len(Q_max))
print('the max-clique is:',Q_max)
print('the run time is:',end-start,'s')


# brock200_1.txt




No = {}
Q = []
Q_max = []
E = np.loadtxt("brock200_1.txt" ,dtype=int) 
V = set(E.flatten()) 
start = time.time()
Rmin,VV = SORT(E,V)
Rmin,No = NUMBER_SORT(E,V,Rmin,No)
for i in range(1,len(Rmin)):
    VV[i] = Rmin[i]
Q_max = NUMBER (E,V, Rmin, VV)
end = time.time()
print('the max-clique cardinality is:',len(Q_max))
print('the max-clique is:',Q_max)
print('the run time is:',end-start,'s')


# brock200_2.txt




No = {}
Q = []
Q_max = []
E = np.loadtxt("brock200_1.txt" ,dtype=int) 
V = set(E.flatten()) 
start = time.time()
Rmin,VV = SORT(E,V)
Rmin,No = NUMBER_SORT(E,V,Rmin,No)
for i in range(1,len(Rmin)):
    VV[i] = Rmin[i]
Q_max = NUMBER (E,V, Rmin, VV)
end = time.time()
print('the max-clique cardinality is:',len(Q_max))
print('the max-clique is:',Q_max)
print('the run time is:',end-start,'s')


# p_hat300_1.txt




No = {}
Q = []
Q_max = []
E = np.loadtxt("p_hat300_1.txt" ,dtype=int) 
V = set(E.flatten()) 
start = time.time()
Rmin,VV = SORT(E,V)
Rmin,No = NUMBER_SORT(E,V,Rmin,No)
for i in range(1,len(Rmin)):
    VV[i] = Rmin[i]
Q_max = NUMBER (E,V, Rmin, VV)
end = time.time()
print('the max-clique cardinality is:',len(Q_max))
print('the max-clique is:',Q_max)
print('the run time is:',end-start,'s')


# p_hat300_3.txt




No = {}
Q = []
Q_max = []
E = np.loadtxt("p_hat300_3.txt" ,dtype=int)
V = set(E.flatten()) 
start = time.time()
Rmin,VV = SORT(E,V)
Rmin,No = NUMBER_SORT(E,V,Rmin,No)
for i in range(1,len(Rmin)):
    VV[i] = Rmin[i]
Q_max = NUMBER (E,V, Rmin, VV)
end = time.time()
print('the max-clique cardinality is:',len(Q_max))
print('the max-clique is:',Q_max)
print('the run time is:',end-start,'s')








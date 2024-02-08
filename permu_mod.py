class permu_modules:
    def __init__(self,mat,sub):
        self.mat=mat
        self.sub=sub
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
def regular_rep(G):
    g=G.gens()
    sub=G.all_subgroups()
    quot=[[] for _ in range(len(sub))]
    bas=[[] for _ in range(len(sub))]
    for j in range(len(sub)):
        quot[j]=G.quotient(sub[j])
        bas[j]=quot[j].list()

    L=[[] for _  in range(len(sub))]
    for j in range(len(sub)):
        L[j]=[[] for _ in range(len(g))]
        for i in range(len(g)):
            L[j][i]=[[] for _ in range(quot[j].order())]
    
            for k in range(quot[j].order()):
                L[j][i][k]=[0 for j in range(quot[j].order())]
                L[j][i][k][bas[j].index(quot[j](g[i])*bas[j][k])]=1

    mat=[[matrix(ZZ,L[j][i]).transpose() for i in range(len(g))] for j in range(len(sub))]
    #ge= [x.order() for x in g ] + [p]
    return permu_modules(mat,sub)
def diag_permu(G):
    p=prime_divisors( G.order())[0]
    F = pAdicField( p, 20, print_mode = "digits" )
    reg=regular_rep(G)
    mats=reg.mat
    mats.remove(mats[0]), mats.remove(mats[len(mats)-1])
    mats=[[matrix(F,x) for x in mats[j]]  for j in range(len(mats))]
    diag=[butler_diagram(G,mats[i]) for i in range(len(mats))]
    return diag

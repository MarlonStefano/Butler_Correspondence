class permu_modules:
    def __init__(self,mat,sub,p):
        self.mat=mat
        self.sub=sub
        self.p=p
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
def regular_rep(G):
    g=G.gens()
    sub=G.all_subgroups()
    p=prime_divisors( G.order())[0]
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
    
            for k in range(len(quot[j].order())):
                L[j][i][k]=[0 for j in range(len(quot[j].order()))]
                L[j][i][k][bas[j].index(quot[j](g[i])*bas[j][k])]=1

    mat=[[matrix(ZZ,L[j][i]).transpose() for i in range(len(g))] for j in range(len(sub))]
    #ge= [x.order() for x in g ] + [p]
    return permu_modules(mat,sub,p)

#G=AbelianGroupGap(regular_rep(G).ge)  # G\simeq H\times C_p           
#mats=regular_rep(H).mat +[identity_matrix(ZZ,(regular_rep(H).p**-1)*G.order())]# incompleto

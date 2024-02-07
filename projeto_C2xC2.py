class Lattices_C2xC2:
    def __init__(self,act_latt, mats, gen_Vi):
        self.act_latt=act_latt
        self.mats=mats
        self.gen_Vi=gen_Vi
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
def Diag_C2xC2(k):
    V=AbelianGroupGap([2]*k)
    sub=V.all_subgroups()
    L2=[[] for _ in range(len(sub)) ]
    for i in range(len(sub)):
        L2[i]=[[sub[i],sub[j]] for j in range(len(sub))]
    L3=[[] for _ in range(len(sub))]
    for i in range(len(sub)):
        L3[i]=[[L2[i][j]+[sub[k]] for j in range(len(sub)) ] for k in range(len(sub))] 
    L4=[[] for _ in range(len(sub))]
    for i in range(len(sub)):
        L4[i]=[[[L3[i][k][j]+[sub[m]] for k in range(len(sub))] for j in range(len(sub))] for m in range(len(sub))]    
    L4=sum(sum(sum(L4, [] ), []), [])

    di_gen=[[] for _ in range(len(L4))]
    gen_che=[[] for _ in range(len(L4))]
    for i in range(len(L4)):
        di_gen[i]=[[x for x in L4[i][j].gens() ] for j in range(4)]
        gen_che[i]=[[di_gen[i][j] for j in range(4) if j!=m] for m in range(4)]
        gen_che[i]=[sum(gen_che[i][j],[]) for j in range(4)]
    diag=[ ]
    for i in range(len(L4)): 
        if all(V.subgroup(gen_che[i][j]).order()==V.order() for j in range(4)):
           diag.append(L4[i])          
    gen_di=[[diag[i][j].gens() for j in range(4)] for i in range(len(diag))]
    gen_Vi=[[[ V(x).exponents() for x in gen_di[i][j]] for j in range(4)] for i in range(len(diag))]
    for j in range(len(gen_Vi)):
        for i in range(4):
            if len(gen_Vi[j][i])==0:
                gen_Vi[j][i]=[vector([0 for j in range(k)])]
    mats=[[] for _ in range(len(gen_Vi))]
    for i in range(len(gen_Vi)):
        mats[i]=[[identity_matrix(IntegerModRing(2),len(gen_Vi[i][j])),identity_matrix(IntegerModRing(2),len(gen_Vi[i][j]))] for j in range(4)]
                
    return diag, gen_Vi, mats
def lattices_C2xC2(k):
    d=Diag_C2xC2(k)
    mats=d[2]
    gen_Vi=d[1]
    m=1
    lamb_info=geral([2,2])
    latt_C2xC2=[[] for _ in range(len(gen_Vi))]
    
    for i in range(len(gen_Vi)):
        latt_C2xC2[i]=[lattice(m,mats[i], d[1][i],lamb_info)]
    act_latt=[latt_C2xC2[i].U_acts for i in range(len(latt_C2xC2))]    
    return Lattices_C2xC2(act_latt, mats, gen_Vi)

                              

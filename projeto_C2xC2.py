class Lattices_C2xC2:
    def __init__(self,act_latt, mats, gen_Vi):
        self.act_latt=act_latt
        self.mats=mats
        self.gen_Vi=gen_Vi
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
def Diag_C2xC2(k):
    W=VectorSpace(IntegerModRing(2),k)
    M=MatrixSpace(IntegerModRing(2),k)
    V=AbelianGroupGap([2]*k)
    aut=V.aut()
    sub=V.all_subgroups()
    L=[[x,y,z,w] for x in sub for y in sub for z in sub for w in sub]

    di_gen=[[] for _ in range(len(L))]
    gen_che=[[] for _ in range(len(L))]
    for i in range(len(L)):
        di_gen[i]=[[x for x in L[i][j].gens() ] for j in range(4)]
        gen_che[i]=[[di_gen[i][j] for j in range(4) if j!=m] for m in range(4)]
        gen_che[i]=[sum(gen_che[i][j],[]) for j in range(4)]
    diag=[ ]
    for i in range(len(L)): 
        if all(V.subgroup(gen_che[i][j]).order()==2**k for j in range(4)):
           diag.append(L[i])
    gen_di=[[diag[i][j].gens() for j in range(4)] for i in range(len(diag))]
    gen_Vi=[[[ V(x).exponents() for x in gen_di[i][j]] for j in range(4)] for i in range(len(diag))]
    for j in range(len(gen_Vi)):
        for i in range(4):
            if len(gen_Vi[j][i])==0:
                gen_Vi[j][i]=[vector([0 for j in range(k)])]
    pro=[x for x in M if x**2==x]
    pro_diag=[[x for x in pro if  all(all(x*W(w) in W.subspace([y for y in gen_Vi[i][j]]) for w in gen_Vi[i][j]) for j in range(4)) ] for i in range(len(gen_Vi))   ]
    idecom=[gen_Vi[i] for i in range(len(gen_Vi)) if len(pro_diag[i])==2]
    idecom_gro=[diag[gen_Vi.index(x)] for x in idecom] 

    cl=[[] for _ in range(len(idecom_gro))]
    for i in range(1,len(idecom_gro)):
        cl[0]=[x for x in idecom_gro if any(all(V.subgroup([f(idecom_gro[0][j].gens()[l]) for l in range(len(idecom_gro[0][j].gens()))]).is_subgroup_of(x[j])  and V.subgroup([f(idecom_gro[0][j].gens()[l]) for l in range(len(idecom_gro[0][j].gens()))]).order()==x[j].order() for j in range(4) ) for f in aut)]
        cl[i]=[x for x in idecom_gro if any(all(V.subgroup([f(idecom_gro[i][j].gens()[l]) for l in range(len(idecom_gro[i][j].gens()))]).is_subgroup_of(x[j])  and V.subgroup([f(idecom_gro[i][j].gens()[l]) for l in range(len(idecom_gro[i][j].gens()))]).order()==x[j].order() for j in range(4) ) for f in aut) and all(x not in cl[k] for k in range(i)) ]
    cl=[x for x in cl if x!=[]] 
    cl=[cl[j][0] for j in range(len(cl))]
    cl_gen_Vi=[[[ V(x).exponents() for x in cl[i][j].gens()] for j in range(4)] for i in range(len(cl))]
    
    mats=[[] for _ in range(len(gen_Vi))]
    for i in range(len(gen_Vi)):
        mats[i]=[[identity_matrix(IntegerModRing(2),len(gen_Vi[i][j])),identity_matrix(IntegerModRing(2),len(gen_Vi[i][j]))] for j in range(4)]

                                         
    return diag, gen_Vi, mats,cl,pro_diag,idecom,cl_gen_Vi
def lattices_C2xC2(k):
    d=Diag_C2xC2(k)
    mats=d[2]
    gen_Vi=d[1]
    m=1
    lamb_info=geral([2,2])
    latt_C2xC2=[lattice(m,mats[i], d[1][i],lamb_info) for i in range(len(d[1]))]
    act_latt=[latt_C2xC2[i].U_acts for i in range(len(latt_C2xC2))]    
    return Lattices_C2xC2(act_latt, mats, gen_Vi)



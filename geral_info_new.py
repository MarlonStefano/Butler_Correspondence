class lambda_i:
    """As informações relacionadas a ação de G em cada Lambda_i
    self.p:               o primo p
    self.gen_group:       os geradores do grupo G
    self.div_ele:         a lista de divisores elementares de G
    self.sub:              a lista de subgrupos de G para o qual o quociente é cílico
    self.gen_quoti:       a lista de geradores dos quocientes cíclicos de G
    self.mat_gen          a lista da lista de matrizes da ação de cada gerador de G nos Lambdas_i
    self.gen:             a lista de geradores de G que também geram os quocientes cíclicos
    """
    def __init__(self,p,gen_group,div_ele,sub,gen_quoti0,mat_gen,gen,G,quoti_order):
        self.p=p              
        self.gen_group=gen_group       
        self.div_ele=div_ele
        self.sub=sub             
        self.gen_quoti=gen_quoti0
        self.mat_lambda_i=mat_gen
        self.gen=gen
        self.G=G
        self.quoti_order=quoti_order
def __repr__(self):
    return f' As informções da ação de {self.G} nos Lambdas_i'                  
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
# função que fornece a matriz da ação de G em Lambda_i através de um gerador específico de um quociente de G que é cíclico. 
#toma como input um primo p e a ordem do quociente.

# função que devolve uma lista de geradores de cada quociente de G que é cíclico e sua matriz de ação em Lambda_i.
#  Toma como input uma lista de divisores elementares de G.

def geral(div_ele):
    
    
    G=AbelianGroupGap(div_ele)
    p=ZZ(prime_divisors(G.order())[0])
    sub=[ x for  x in subgroups_with_cyclic_quotient( G ) if x.order()<G.order()]
    gen_group=gens(G)
    gen_quoti0=[[] for _ in range(len(sub))]
    for j in range(len(sub)):
        gen_quoti0[j]=[G.quotient(sub[j]).lift(x) for x in gens(G.quotient(sub[j])) if x.order()==G.quotient(sub[j]).order()]
        gen_quoti0[j]=gen_quoti0[j][0]
    
    P=[[] for _ in range(len(gens(G)))]
    for j in range(len(gens(G))):
        P[j]=[[x for x in range(G.quotient(sub[i]).order()) if G.subgroup([gens(G)[j]*gen_quoti0[i]**x]).is_subgroup_of(G.subgroup(sub[i]))==true] for i in range(len(gen_quoti0))]
        P[j]=[P[j][k][0] for k in range(len(P[j]))]
    Pw=[[] for _ in range(len(gens(G)))]
    for j in range(len(gens(G))):
        y=var('y')
        Pw[j]=[[x for x in solve_mod(y==-P[j][k],G.quotient(sub[k]).order() )[0]] for k in range(len(P[j]))] 
        Pw[j]=[Pw[j][k][0] for k in range(len(Pw[j]))]
        mat_gen=[[] for _  in range(len(gens(G)))]
    for j in range(len(gens(G))):
        mat_gen[j]=[matrix(ZZ,[[1]])]+[ mat_lam(p,G.quotient(sub[i]).order())**Pw[j][i] for i in range(len(Pw[j]))] 
    gen=[[] for _ in range(len(sub))]
    for j in range(len(sub)):
        gen[j]=[x for x in gens(G) if G.subgroup([x*gen_quoti0[j]**(G.quotient(sub[j]).order()-1)]).is_subgroup_of(G.subgroup(sub[j]))==true]
        gen[j]=gen[j][0]
    quoti_order=[G.quotient(sub[k]).order() for k in range(len(sub))]            
    return lambda_i(p,gen_group,div_ele,sub,gen_quoti0,mat_gen,gen,G,quoti_order)

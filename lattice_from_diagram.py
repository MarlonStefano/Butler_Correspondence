## Código escrito por Marlon
## Modificado por Csaba
## Modificado por Marlon novamente para trabalhar com um p-grupo finito abeliano arbitrário

class Reticulado_ass_diagrama:
    """The data structure to store information related to a GZ_p-lattice."""
    def __init__( self, G,U_acts,mats):
    
        self.group = G
        self.U_acts=U_acts
        self.mats = mats 
        self.p = prime_divisors( G.order())[0]
        self.padic_ring = mats[0][0][0].parent()
         
# mat is matrix
# gen_U is generating set of Z-module. 
# U invariant under mat
# returns the matrix of the action of mat on U

def act_U(mat,gen_U):
    l=[mat*gen_U[j] for j in range(len(gen_U))]
    m=transpose(matrix(ZZ,gen_U))
    coor=[[] for _ in range(len(l))]
    for j in range(len(l)):
        coor[j]=m.solve_right(vector(ZZ,l[j]))
    mact=transpose(matrix(ZZ,coor))
    return mact

# p é primo
# k >= 1
# gen1 é lista de vetores com entradas em Z/(p^k) de comprimento l
# gen2 é lista de vetores com entradas em Z/(p^k) de comprimento l
# os elementos de gen1 pertencem ao módulo gerado por gen2
# mat é matriz quadrada m x m onde m == len( gen2 ) com entradas em Z/(p^k)
#
# output: devolve uma lista das imagens dos vetores em gen1 pelo morfismo V -> V determinado pela matriz mat 

def mudança(p,k,gen1,gen2,mat):# gen1 é lista de elementos de um módulo que desejamos conhecer a imagem por mat e gen2 os geradores desse módulo que ja conhecemos a imagem por mat
    l=len(vector(gen1[0]))
    F=ZZ**l
    R=IntegerModRing(p**k)
    mgen=transpose(matrix(R,[x for x in gen2]))

    m=[[] for _ in range(len(gen1))]
    for j in range(len(gen1)):
        m[j]=mgen.solve_right(vector(R,gen1[j])) 

    coorde=[[] for _ in range(len(gen1))]
    for j in range(len(gen1)):
        coorde[j]=mat*m[j]

    im=[[] for _ in range(len(gen1))]
    for j in range(len(gen1)):
        im[j]=sum(coorde[j][i]*F(gen2[i]) for i in range(len(gen2)))
    return im
# função que fornece a matriz da ação de G em Lambda_i através de um gerador específico de um quociente de G que é cíclico. 
#toma como input um primo p e a ordem do quociente.
def mat_lam(p,indi):
     l=p**(-1)*indi*(p-1)
     L=[[] for _  in range(l)]
     for j in range(l):
        L[j]=[0 for i in range(l)]
     for j in range(l-1):
        L[j][j+1]=1
     for j in range(p-1):
        L[l-1][p**(-1)*indi*j]=-1
     return transpose(matrix(ZZ,L))

def morfis(geradores, g): 
    F=ZZ**len(geradores[0])
    m=[[] for _ in range(len(g.rows()[0]))]
    for j in range(len(g.rows()[0])):
        m[j]=sum(g[i][j]*F(geradores[i]) for i in range(len(g.rows()[0])))

    return m


# create the lattice from the butler diagram

def lattice_from_diagram(diag,lamb):

    #G = diag.group
    #ids = diag.idempotents
    #k = -depth_list( [ x.coefficients() for x in ids ])
    #r = len( diag.Vi )
    #mats = diag.action_Vi
    #V_i_gens = diag.Vi
    #mat_lambda_i = [ Lambda_i_basis( x ) for x in ids ]


    # primeiro recolhemos as informações básicas sobre os Lambda_i e a ação de G em cada um.

    
    # Aqui calculamos JVi e os quocientes Vi/JVi
    k=diag.k
    mats=diag.action_Vi
    V_i_gens=[diag.Vi[j].rows() for j in range(len(diag.Vi))]
    G=lamb.G
    r=len(lamb.sub)+1
    p=lamb.p
    R=Zp(p)

    l=len(vector(V_i_gens[0][0]))
    
    F0=ZZ**l
    #return F0

    F1=F0.submodule(p**k*gens(F0)[j] for j in range(l))

    F2=F0/F1

    Vi=[[] for _ in range(r)]

    for j in range(r):
        Vi[j]=F2.submodule([F2(x) for x in V_i_gens[j] ]) 
    
    JV_gens=[[] for _ in range(r)]

    for j in range(r): #versão anterior
        JV_gens[j]=[[F2(morfis(V_i_gens[j],mats[j][i])[m])-F2(V_i_gens[j][m]) for m in range(len(V_i_gens[j]))]  for i in range(len(lamb.gen_group))]+[[p*F2(x) for x in gens(Vi[j])]] # a ação de g deve ser dada por um morfismo de Vj
   
     

    JV=[[] for _ in range(r)]

    for j in range(r):
        JV[j]=F2.submodule([F2(x) for x in sum( JV_gens[j], [] )])
     
    W=[[] for _ in range(r)]         # lista dos quocientes V_i/JV_i

    D=[]               # listas das dimensões

    for j in range(r): 
        W[j]= Vi[j].V()/JV[j].V()
        d=len(gens(W[j]))
        D.append(d)
    
    s=  D[0]+sum( D[j+1]*(p-1)*p**(-1)*lamb.quoti_order[j] for j in range(len(lamb.sub)))    

    W_gens=[[] for _ in range(r)]

    for j in range(r):
        W_gens[j]=[x for x in gens(W[j])]

    W_lift=[[] for _ in range(r)]

    for j in range(r):
        W_lift[j]=[W[j](x).lift() for x in W_gens[j]]
 
    V_pro=[[] for _ in range(r)] # lista das listas dos levantamentos dos geradores de V_i/JV_i
 
    for j in range(r):
        V_pro[j]=[F2(x) for x in W_lift[j]]
     
    Vi_f_im=[[] for _ in range(len(lamb.sub))]
 
    for j in range(len(lamb.sub)):
        Vi_f_im[j]=[[mudança(p,k, [V_pro[j+1][m]], V_i_gens[j+1], mats[j+1][lamb.gen_group.index(lamb.gen[j])]**i)  for i in range(p**(-1)*lamb.quoti_order[j]*(p-1))] for m in range(D[j+1])]
    Vi_f_im=sum(sum(sum(Vi_f_im, [] ), []), [])   
    V0_f_im=[x for x in V_pro[0]]
 
    V_f_im=V0_f_im  + Vi_f_im
              
    F3=ZZ**s
    F4=F3.submodule([p**k*gens(F3)[j] for j in range(s)])
    F5=F3/F4
    f=F5.hom([F2(x) for x in V_f_im])
    U0=f.kernel()
    U0_gens=[F5(x) for x in gens(U0)]
    U_gens=[F5(x).lift() for x in U0_gens]+[p**k*y for y in gens(F3)]
    U_mat=matrix(ZZ,U_gens).echelon_form(include_zero_rows=False)
    U=U_mat.rows()
    
    F_act=[[] for _ in range(len(lamb.gen_group))] # retorna erro
    for j in range(len(lamb.gen_group)):
        blocks = sum( [D[x]*[lamb.mat_lambda_i[j][x]] for x in range(len(D))] , [] )
        F_act[j] = block_diagonal_matrix( blocks ) #mat_act(D,mat_lambda_i[j],r)
        
    U_acts=[[] for _ in range(len(lamb.gen_group))]
    for j in range(len(lamb.gen_group)):
        U_acts[j]=act_U(F_act[j],U)

     
    return Reticulado_ass_diagrama(G, U_acts,mats)
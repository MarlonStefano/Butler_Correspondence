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

# create the lattice from the butler diagram given by a lattice

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
    k=-diag.depth
    p=lamb.p
    mats= [[matrix(IntegerModRing(p**k), diag.action_Vi[j][i].rows()) for i in range(len(diag.action_Vi[j]))] for j in range(len(diag.action_Vi))]
    V_i_gens=[diag.Vi[j].rows() for j in range(len(diag.Vi))]
    G=lamb.G
    r=len(lamb.sub)+1
    
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
def leading_position( vec ):
    return vec.support()[0]
    
def depth( el ):
    """The depth of a vector with entries in a p-adic field or integer ring. 
       Returns k where p^k is the minimum valuation among the entries in el."""
    return min( valuation(c) for c in el )

def depth_list( ellist ):
    """Returns the minimum among the depths of ellinst."""
    return min( depth( el ) for el in ellist )

def depth_matrix( mat ):
    """mat is a matrix with entries in a p-adic field or integer ring. Returns the 
    minimum among the valuations of the entries of mat."""
    r, c = mat.nrows(), mat.ncols()
    return depth( vector([ mat[i,j] for i in range(r) for j in range( c )])) 

def depth_matrix_list( mlist ):
    """return the minimum among the depths of the matrices in mlist"""
    return min( depth_matrix( mat ) for mat in mlist )


# compute the complement of U in W
def complement( W, U ):
    
    # check if subspace 
    assert U.is_subspace( W )

    # calculate basis for U
    vectsU = [ u for u in U.basis()] 

    # this will hold the generating set of the complement   
    vects = []
    for v in W.basis():
        if not (v in U):        
            vects.append( v )
            vectsU.append( v )
            U = W.subspace( vectsU )

    return W.span( vects )


def image_vector_under_action( V, mat, vec ):

    # get the coefficients of vec in the linear combination 
    # of the rows of mat
    return V.solve_left( vec )*mat*V 
#funções para lattice.py 
def ac_U(mat,gen_U,R0):
    l=[mat*gen_U[j] for j in range(len(gen_U))]
    m=transpose(matrix(R0,gen_U))
    coor=[[] for _ in range(len(l))]
    for j in range(len(l)):
        coor[j]=m.solve_right(vector(R0,l[j]))
    mact=transpose(matrix(R0,coor))
    return mact

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
def right_regular_action( g ):

    g_elts = g.list()
    g_gens = g.gens()
    perms = []
    sym_gp = SymmetricGroup( g.order()) 
    for g in g_gens:
        perm = sym_gp( [ g_elts.index( x*g ) + 1 for x in g_elts ])
        perms.append( perm )
    
    return sym_gp.subgroup( perms, canonicalize = False )

# funções para idempotents.py
def subgroups_with_cyclic_quotient( g ):
    gp = right_regular_action( g )

    subs_gp = [ x for x in gp.subgroups() if gp.quotient(x).is_cyclic() ]
    no_gens = len( g.gens())
    gen_dict = { str(gp.gens()[k]): g.gens()[k] for k in range( no_gens )}
    #return gen_dict
    
    subs = []
    for h in subs_gp:
        gens_h_perm = h.gens()
        gens_h = []
        for x in gens_h_perm: 
            if x == x**0:
                continue
            w = x.word_problem( gp.gens(), as_list = True, display = False )
            gens_h.append( prod( gen_dict[x[0]]**x[1] for x in w ))
        subs.append( g.subgroup( gens_h ))

    return subs         

    
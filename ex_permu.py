
G = AbelianGroup( div_ele )
R = pAdicField( ZZ(prime_divisors( G.order())[0]) , 3, print_mode = "digits" )
def act_mat_g(G):
    L=[[] for _  in range(G.gens()[0].order())]
    for j in range(G.gens()[0].order()):
        
        L[j]=[0 for i in range(G.gens()[0].order())]
    L[G.gens()[0].order()-1][0]=1
    for j in range(G.gens()[0].order()-1):
        L[j][j+1]=1

        
    return transpose(matrix(R,L))
mats=[act_mat_g(G)]+ [identity_matrix(R,G.gens()[0].order()) for i in range(len(G.gens())-1)]
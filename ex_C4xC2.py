F = pAdicField( 2, 20, print_mode = "digits" )
G = AbelianGroupGap( [4,2] )
mats =[matrix(F,[[1,0,0],[-1,-1,-2],[0,1,1] ]), matrix(F,[[1,0,0],[0,1,0],[0,0,1]])]
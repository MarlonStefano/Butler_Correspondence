def image( el, mats ):
    """el is an element of a group algebra FG and mats is a list of square matrices matrices over F. 
    Substitutes the elements in mats in the place of the generators of G.
    
    Can be used to calculate a homomorphism FG -> M(n,F) for fields where the standard homomorphism 
    machinery of Sage fails."""

    #initialize with the zero matrix
    im = 0*mats[0]
    R_mats = mats[0][0,0].parent()
    
    #take the dictionary with the terms and their coefficients
    terms = el.monomial_coefficients()
    
    for mon, coeff in terms.items():
        exps = mon.exponents()
        pr = prod( mats[m]**exps[m] for m in range(len(exps)))
        im += R_mats(coeff)*pr

    return im 


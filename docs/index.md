---
layout: default
title: Butler Correspondence
---
# Tutorial for using the implementation of the Butler's correspondence
This repository contains an implementation of the Butler's correspondence, which is a bijection between certain combinatorial objects. The code is written in Python and can be used to explore properties of lattices for $\mathbb{Z}_pG$-modules where $G$ is a finite abelian $p$-group and $\Z_p$  is the ring of $p$-adic integers. 
## Theoretical background
Let  $G$ be a finite abelian  $p$-group and let  $\mathbb{Q}_p$ be the field of  $p$-adic numbers. Denote by $e_0, e_1,\ldots, e_r$ the primitive orthogonal idempotents of the semisimple algebra $\mathbb{Q}_p G$. In the article entitled "On the classification of local integral representations of finite abelian p-groups", Butler provides a method to associate most $\mathbb{Z}_G$-lattices to diagrams of the form:
$$
(V; V_0, \ldots, V_r),
$$
where $ V $ is a finitely generated $ \mathbb{Z}_p G $-module with $ p^{-1}|G| V = 0 $, and each $ V_i $ is a $ \mathbb{Z}_p G$-submodule of $ V $, subject to the following conditions:
$ V = \sum_{j\neq i} V_j $ for all $ i = 0, 1, \ldots, r $, and 
if $ e_i $ is a primitive idempotent of $ \mathbb{Q}_p G $ and $ K_i $ is the kernel of the natural ring homomorphism $ \mathbb{Z}_p G \to (\mathbb{Z}_p G)e_i $, then $ K_i V_i = 0 $.


There is a  natural notion of morphism between diagrams, allowing us to establish a category of such objects.

A $ \mathbb{Z}_p G $-lattice is said to be reduced if it does not have a direct summand isomorphic to any of the modules $ (\mathbb{Z}_p G)e_i $ for $ i = 0, \ldots, r $, or to $ \mathbb{Z}_p G $ itself.
Butler's correspondence is defined via a functor that associates a reduced $ \mathbb{Z}_p G $-lattice $ U $ to a diagram $ \Delta(U) = (V; V_0, \ldots, V_r) $. Conversely, each diagram corresponds to a reduced lattice, though this process is not functorial. However, two diagrams are isomorphic if and only if they correspond to isomorphic reduced lattices.

We will briefly describe the process: given a reduced lattice $ U $, we treat $ U $ as its image in $ \mathbb{Q}_p \otimes_{\mathbb{Z}_p} U $. For each idempotent $ e_i $, we associate the lattice $ e_i U \subseteq \mathbb{Q}_p \otimes_{\mathbb{Z}_p} U $. Since $ \sum_{i=0}^r e_i = 1 $, we have $U\leqslant\sum_{i=0}^r e_i U$. In fact, we have:   
$$
U \subseteq \bigoplus_{i=0}^r e_i U \leqslant \mathbb{Q}_p \otimes_{\mathbb{Z}_p} U
.$$
 Then we define the quotients:
$$
V_i = \frac{e_i U + U}{U}, \quad i = 0, \ldots, r, \quad \text{and} \quad V = \sum_{i=0}^r V_i.
$$
Thus, $ \Delta(U) = (V; V_0, \ldots, V_r) $ satisfies the conditions of a Butler diagram. The strength of this correspondence lies in the fact that if $ W \cong U $, then $ \Delta(W) \cong \Delta(U) $.
In the reverse direction, given a Butler diagram $ D = (V; V_0, V_1, \ldots, V_r) $, we construct a reduced $ \mathbb{Z}_p G $-lattice $ L(D) $ as follows. Let $ J $ be the Jacobson radical of $ \mathbb{Z}_p G $. If $ d_i = \dim_{\mathbb{F}_p}(V_i / J V_i) $, for $ i = 0, \ldots, r $, we take the lattice:


$$
F = \bigoplus_{i=0}^r (\Lambda_i)^{d_i}, \quad \text{where } \Lambda_i = (\mathbb{Z}_p G)e_i.
$$
Next, we define a surjective homomorphism $ \varphi : F \to V $, induced by surjective maps $ \varphi_i : (\Lambda_i)^{d_i} \to V_i $. Then we define:


$$
L(D) = \ker(\varphi),
$$
which is a reduced lattice satisfying $ \Delta(L(D)) \cong D $. Moreover, if $ D' \cong D $, then $ L(D') \cong L(D) $. Thus, for every diagram $ D $, there exists a reduced $ \mathbb{Z}_p G $-lattice $ L(D) $ such that $ \Delta(L(D)) \cong D $.
## Installation
To use the implementation, you need to have SageMath installed on your computer. You can download it from the official website: https://www.sagemath.org/download.html. Once you have SageMath installed, you can clone this repository and navigate to the directory containing the code.
## Usage
You can run the code in SageMath by executing the following command in the terminal:

    sage: load(main.py) 

This will load the main.py file, which contains the implementation of the Butler's correspondence. Now you can use the functions defined in the code to explore the diagrams and the correspondents lattices for a fixed finite abelian $p$-group $G$. First, you need to define the group $G$ and get the information about the primitive idempotents  of the group algebra $\Z_pG$. First, we give a list of elementary divisors of $G$, for example, if G is isomorphic to $\Z/2\Z \times \Z/2\Z$, you can define it as the list [2, 2]. So the general information about the group $G$ can be obtained by running the following command:

    sage: div_ele=[2, 2]
    sage: ger=geral(div_ele)

This will give you the general information about the group $G$, including the number of primitive idempotents  and the structure of the modules $(\Z_pG)e_i$, where $e_i$ are the primitive idempotents of $\Z_pG$. So we define $G$ as follows:
   
    sage:G=ger.G

Now you can explore the diagrams and the correspondents lattices for the group $G$. For example, every primitive idempotent $e_i$ corresponds to a subgroup of G such that the quotient of $G$ by this subgroup is cyclic. You can get the list of these subgroups by running the following command:

    sage: sub=ger.sub

This will give you the list of subgroups of $G$ that correspond to the nontrivial primitive idempotents of $\Z_pG$. So every component $V_i$ of a diagram corresponds to a subgroup of $G$ with cyclic quotient and the order that the subgroups are listed in the list sub is the same as the order of the components in the diagrams. For example, for the group $G=\Z/2\Z \times \Z/2\Z$ we can get its generators by running the following command:

    sage: gen=G.gens()

So the list of subgroups that correspond to the non trivial primitive idempotents of $\Z_pG$ is given by:

    sage: sub 
    [Subgroup of Abelian group with gap, generator orders (2, 2) generated by (f2,),
    Subgroup of Abelian group with gap, generator orders (2, 2) generated by (f1,),
    Subgroup of Abelian group with gap, generator orders (2, 2) generated by (f1*f2,)]
 
 So a diagram for the group $G=\Z/2\Z \times \Z/2\Z$ is given by the components $V_0,V_{f2}, V_{f1}, V_{f1*f2}$, where $V_0$ corresponds to the trivial idempotent. The structure of every component $V_i$ is given by the structure of the module $(\Z_pG)e_i$, and the structure of the module $(\Z_pG)e_i$ is given by the matrices of the action of the generators of G on the module $(\Z_pG)e_i$. These matrices can be obtained by running the following command: 
 
    sage: mat=ger.mat_lambda_i
    sage:  mat
    [[[1], [-1], [1], [-1]], [[1], [1], [-1], [-1]]]

This will give you the matrices of the action of the generators of G on the modules (Z_pG)e_i, where the first list corresponds to the action of the first generator on the respective modules $(Z_pG)e_0, (Z_pG)e_{f2}, (Z_pG)e_{f1}, (Z_pG)e_{f1*f2}$ and the second list corresponds to the action of the second generator on these modules. In this case we observe that $G$ must act trivially on the modules $V_i$: one generator acts trivially and the other acts as -1, but since each $V_i$ is a $\mathbb{Z}/2\mathbb{Z}$-module the action by -1 is trivial, so the action of $G$ on each $V_i$ is trivial.  This gave us the structure of the diagrams for the group $G=\Z/2\Z \times \Z/2\Z$.

## Diagrams from lattices
 Now we can proceed to calculate a diagram from a lattice. For this we need to define a lattice, and this is done by defining the matrices of the action of the generators of $G$ on the lattice. For example, the following matrices define a lattice $U$ for the group $G=\Z/2\Z \times \Z/2\Z$:


    sage: F = pAdicField( 2, 20, print_mode = "digits" ) 
    sage: mats=[matrix(F,[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 1, 1, 0], [0, -1, 0, 0, -1]]), matrix(F,[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [-1, 0, -1, -1, 0], [0, -1, -1, 0, -1]])]
 
 As we already defined the group $G$, we can calculate the diagram that corresponds to the 
 lattice $U$ by running the following command:
    
     sage: diag=butler_diagram(G,mats)

This will give you a function that provides information about the components of the diagram, for example the following command will give you the list of matrices such that the i-th matrix is the matrix whose rows are the generators of the component $V_i$ of the diagram, where the order of the components is the same as the order of the subgroups in the list sub. So the first matrix corresponds to the generators of the component $V_0$, the second matrix corresponds to the generators of the component $V_{f2}$, the third matrix corresponds to the generators of the component $V_{f1}$ and the fourth matrix corresponds to the generators of the component $V_{f1*f2}$. So we can get the generators of the components of the diagram by running the following command: 


    sage: diag.Vi 
    [
     [0 0 0 1 0]
     [0 0 0 0 1], [0 0 0 1 1], [0 0 0 1 0], [0 0 0 0 1]].

The module $V$ which is the sum of the components $V_i$ of the diagram is given by the following command:

    sage: diag.V
    [0 0 0 1 0]
    [0 0 0 0 1]

This is a matrix whose rows are the generators of the module $V$. We can get the matrices of the action of the generators of $G$ on the module $V$ by running the following command:

    sage: diag.action_V
    [
     [1 0]  [1 0]
     [0 1], [0 1]
    ]

This is a list of matrices where the first matrix corresponds to the action of the first generator of $G$ on the module $V$ and the second matrix corresponds to the action of the second generator of $G$ on the module $V$. What we mean by the matrices of the action is that these matrices are the matrices of the action of $G$ relative to the generators of $V$ that appear in the matrix diag.V. The matrices of the action of $G$ on the components $V_i$ of the diagram can be obtained by running the following command:

    sage: diag.action_Vi
    [[
    [1 0]  [1 0]
    [0 1], [0 1]
    ], [[1], [1]], [[1], [1]], [[1], [1]]] 

This is a list of lists of matrices, where the $i$-$th$ list corresponds to a list of matrices that give the action of the generators of $G$ on the component $V_i$ of the diagram, The matrices are the matrices of the action of $G$ relative to the generators of $V_i$ that appear in the $i$-$th$ matrix of diag.Vi. So for example the first list corresponds to the component $V_0$ and the first matrix in this list corresponds to the action of the first generator of $G$ on $V_0$.

### Permutation lattices
We can calculate the correspondents diagrams for indecomposables permutations lattices over a finite abelian $p$-group $G$. First we need to define $G$ as we did before. This time we will work with the group $G=\Z/4\Z \times \Z/2\Z$ for make explicit some differences with the previous example. So we can define $G$ as follows:

    sage: div_ele=[4, 2]
    sage: ger=geral(div_ele)
    sage:G=ger.G

Now we can get the list of diagrams that correspond to the permutation lattices over $G$ by running the following command:

    sage: diag_permu=diag_permu(G)

Every indecomposable permutation module corresponds to a subgroup of $G$, and the order of the subgroups of $G$ on G.all_subgroups() is the same as the order of the diagrams on the list diag_permu(G) except the first subgroup and the last subgroup of G.all_subgroups() which correspond to the trivial subgroup and the whole group $G$ respectively. So for example we have the following list of subgroups of G:

    sage: gen=G.gens()
    sage gen
    (f1, f3)
    sage: G.all_subgroups()
    [Subgroup of Abelian group with gap, generator orders (4, 2) generated by (),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f3,),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f2,),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f2*f3,),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f2, f3),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f1, f2),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f1*f3, f2),
    Subgroup of Abelian group with gap, generator orders (4, 2) generated by (f2, f3, f1)]

Hence, the first diagram in the list diag_permu(G) corresponds to the diagram of the permutation module $\Z_2[G/<f3>]$ which is:

    sage: diag_permu[0].Vi
    [
               [2 0 2 0]
    [1 1 1 1], [0 2 0 2], [], [1 3 1 3], [], []] 
    sage: diag_permu[0].action_Vi
    [[[1], [1]],
    [
    [0 1]  [1 0]
    [1 0], [0 1]],
    [[], []], [[3], [1]], [[], []], [[], []]]
    sage: diag_permu[0].V
    [1 1 1 1]
    [0 2 0 2]
    sage: diag_permu[0].action_V
    [
    [1 2]  [1 0]
    [0 3], [0 1]
    ]

To know the smallest quotient of $\Z_p$ such that $V$ is a module over this quotient we can run the following command:

    sage: diag_permu[0].residue_class_ring
    Ring of integers modulo 4

   We also can recover the initial lattice from the last diagram by running the following command:

    sage: U=lattice_from_diagram(diag_permu[0], ger)

The matrices of the action of the generators of $G$ on the lattice $U$ can be obtained by running the following command:

    sage: U.U_acts
    [
    [ 1  0  0  0]  [1 0 0 0]
    [-1 -1 -2  0]  [0 1 0 0]
    [ 0  1  1  0]  [0 0 1 0]
    [ 0  0  1 -1], [0 0 0 1]
    ]

The order of the matrices in the list U.U_acts is the same as the order of the generators of $G$. The function lattice_from_diagram() can be used to recover the initial lattice from any diagram that was calculated by the function butler_diagram() or diag_permu(). There also exists a function that allows recovering the diagram from the lattice that was calculated by lattice_from_diagram(); this function is used as follows:

    sage: diag_from= butler_diagram_from_lattice(U)
    sage: diag_from.Vi
    [
               [0 0 2 0]
    [0 0 2 3], [0 0 0 2], [], [0 0 0 1], [], []]

To finish this section, we recall that the information about the structure of the modules $(\Z_pG)e_i$ can be obtained through the function geral() that we used at the beginning to get the information about the group $G$ as in the first example.
## Lattices from diagrams
We can also calculate a lattice from a diagram. Diagrams are given by the components $V_i$ and the action of $G$ on these components. So we can define a diagram by giving the components $V_i$ and these components are given by their generators, so we can define a diagram by giving a list of list where the $i$-$th$ list corresponds to the generators of the component $V_i$. We also need to give the matrices of the action of the generators of $G$ on the components $V_i$, we can give this information by giving a list of lists of matrices, where the $i$-$th$ list corresponds to the list of matrices that give the action of the generators of $G$ on the component $V_i$, and the $j$-$th$ matrix in this list corresponds to the action of the $j$-$th$ generator of $G$ on the component $V_i$. So for example we can define a diagram for the group $G=\Z/2\Z \times \Z/2\Z$ as follows:

    sage: k=1 # the exponent k such |G|/p=p^k
    sage: div_ele=[2, 2]
    sage: ger=geral(div_ele)
    sage: V_i_gens=[[[1,0],[0,1]],[[1,1]],[[1,0]],[[0,1]]] # the generators of the components V_i of the diagram
    sage: mats=[[matrix(IntegerModRing(2),[[1,0],[0,1]]),matrix(IntegerModRing(2),[[1,0],[0,1]])],[matrix(IntegerModRing(2),[[-1]]),matrix(IntegerModRing(2),[[1]])],[matrix(IntegerModRing(2),[[1]]),matrix(IntegerModRing(2),[[-1]])],[matrix(IntegerModRing(2),[[-1]]),matrix(IntegerModRing(2),[[-1]])]] # the matrices of the action of the generators of G on the components V_i of the diagram

Now we can calculate the lattice that corresponds to this diagram by running the following command:

    sage: U=lattice(k, mats, V_i_gens, ger)
    sage: U.U_acts
    [
    [ 1  0  0  0  0]  [ 1  0  0  0  0]
    [ 0  1  0  0  0]  [ 0  1  0  0  0]
    [ 0  0 -1  0  0]  [ 0  0  1  0  0]
    [ 0  0  1  1  0]  [-1  0 -1 -1  0]
    [ 0 -1  0  0 -1], [ 0 -1 -1  0 -1]
    ]

U.U_acts is a list of matrices where the first matrix corresponds to the action of the first generator of $G$ on the lattice $U$ and the second matrix corresponds to the action of the second generator of $G$ on the lattice $U$. In general, the order of the matrices in the list U.U_acts is the same as the order of the generators of $G$. As earlier, we can recover the diagram from the lattice $U$ by running the following command:

    sage: diag_from= butler_diagram_from_lattice(U)
    sage: diag_from.Vi
    [
    [0 0 0 1 0]
    [0 0 0 0 1], [0 0 0 1 1], [0 0 0 1 0], [0 0 0 0 1]
    ]
    sage: diag_from.action_Vi
    [[
    [1 0]  [1 0]
    [0 1], [0 1]
    ], [[1], [1]], [[1], [1]], [[1], [1]]]

Recall that to define a diagram we need first to get the information about the group $G$ as the list of subgroups of $G$ that correspond to the non trivial primitive idempotents of $\Z_pG$ and the matrices of the action of the generators of $G$ on the modules $(\Z_pG)e_i$, and this information can be obtained by running the function  geral() as we did at the beginning of this tutorial.   

# Dominoes
---
### Questions
 *  Use valarray for operations  
 .. valarray.apply doesn't take capturing lambdas -> workaround?
 .. is it worth it? performance compared to vector
 *  Have GslQuad object as static member or initilize it newly each time a
    corresponding function is called?  
 .. How to implement a static solution? How to pass a member method without
 too much wrapping?
 *  What about the destructor of `gsl_integration_workspace`
 [here](https://stackoverflow.com/a/24151084/9133910)? Why is it better than
 passing `gsl_workspace_free` as a functor?

---
### API
##### Intrinsic speed
 *  Input  
 ... 1d-vector λ (spacing)  
 ... scalar N (pieces to be considered)  
 ... μ (friction)  
 ... `Domino` struct  

 *  Output  
 ... 2d-vector holding φ'\* and V\*  

##### Speed at x
 *  Input  
 ... scalar D (total number of dominoes)  
 ... scalar λ (spacing)
 ... scalar N (pieces to be considered)  
 ... μ (friction)  
 ... `Domino` struct  

 *  Output  
 ... 3d-vector holding x, φ' and V  

---
### Implementation
##### type Domino (struct)
 *  length L
 *  width h

##### DominoChain class
 *  constructor  
 ... m_integrator = Integrator  
 ... m\_N = N  
 ... m\_D = D  
 ... m\_μ = μ  
 ... m\_L = domino.length  
 ... m\_h = domino.width  
 ... m\_λ = λ  
 ... m\_λs = λ_vec  
 ... φ = arctan(h/L)  
 ... ω = √(3*g*cos(φ)/(2*L))  
 ... ξ = L*cos(ψ)  
 ... ψ = ψ(λ)  
 ... θ^ = θ^(λ)  
 ... R = R(λ)  
 ... η = η(λ)  

 *  public  
 ... makeVecFromVec(func, λs, params, result)  
 ... makeVecFromNum(func, D, params, result)  
 ... intrinsicAngular(λ)  
 ... intrinsicTransversal(λ, intrinsicAngularVal)  
 ... angularAtX(i, startVal)  
 ... transversalAtX(i, startAngularVals)  
 ... thetaDot(θ, λ, φ')

 *  private
 // f(λ)
 ... ψ = arcsin(λ/L)  
 ... θ^ = arccos(h/(h+λ))  
 ... R = 1+(ξ+μ*λ)/(ξ-μ*h)  
 ... η = (λ+h)/L
 // 
 ... [P/K](θ, φ')
 ... k(θ, λ)
 ... [θ'_i/θ'_(i+1)](θ, λ)
 ... θ_next(θ, λ)

##### Integrator class
// wrapper for gsl integration  
// see: [scicomp link](https://scicomp.stackexchange.com/a/27248)  
// and: [so link](https://stackoverflow.com/a/24151084/9133910)

---
### Program flow
 *  calculate φ'(x) and V(x) seperatly even though V requieres φ', therewidth
    leading to calculating k(θ=0) twice  
    Reason: calculating these quantities simultaniously would requiere passing
    a result vector reference for φ' to V and accessing it at each integration
    step. The possibly gained performance I do not consider worthy of its price
    of a large increase in complexity (which I cannot estimate enough at this
    point) and a decrease of readability of the code  

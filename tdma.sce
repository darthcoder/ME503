//First we define an array of deltax values.
del = [0.2 0.1 0.05 0.025]

//Now N is 1/deltax
for i = 1:4            // 4 cases of deltax.
    N = 1/del(i)
    
    //Now we define the equations in tridiagonal form 
    a(1) = 0
    b(1) = -2
    c(1) = 1
    d(1) = 3*del(i)^2 - 273
    a(N) = 1
    b(N) = -2
    c(N) = 0
    d(N) = 3*del(i)^2 - 318
    for k = 2:N-1
        a(k) = 1
        b(k) = -2
        c(k) = 1
        d(k) = 3*del(i)^2
    end
    
    //Now we solve
    
    //Define gam(1) and bet(1)
    gam(1) = d(1)/b(1)
    bet(1) = b(1)
    
    //Define rest of the gamma(j) and beta(j)
    for j = 2:N
        bet(j) = b(j) - (a(j)*c(j-1))/bet(j-1)
        gam(j) = (d(j) - a(j)*gam(j-1))/(b(j)-a(j)*c(j-1)/bet(j-1))
    end
    
    //Now to find T
    T(N) = gam(N)
    for m = N-1:-1:1
        T(m) = gam(m) - c(m)*T(m+1)/bet(m)
    end
    
    //storing at a common place
    for n = 1:N
        temp(i,n) = T(n)
    end
    
end
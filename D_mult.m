function D = D_mult(sys,x,u)
    x_mult = reshape(sys.Ap,[sys.n*sys.p,sys.n])*x;
    D= reshape(x_mult,[sys.n,sys.p]);
end
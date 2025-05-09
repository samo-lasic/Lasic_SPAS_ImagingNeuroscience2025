function f=DwSpherical(w,R,D0,alfa,n)

Ks=BesselKernelsSphere(n);

B=2*(R./Ks).^2./(Ks.^2-2);
a=(Ks/R).^2;

for i=1:length(w)
    f(i)=D0*alfa+D0*(1-alfa)*sum(a.*B*w(i)^2./(a.^2*D0^2+w(i)^2));
end


function f=BesselKernelsSphere(n)

for l=1:n
    f(l)=fzero(inline('x*besselj(1/2,x)-2*besselj(3/2,x)'),2+(l-1)*pi);
end

f=f';
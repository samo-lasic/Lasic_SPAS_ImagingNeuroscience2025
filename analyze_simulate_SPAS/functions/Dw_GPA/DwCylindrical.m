function f=DwCylindrical(w,R,D0,alfa,n)

Kc=BesselKernelsCylinder(n);


B=2*(R./Kc).^2./(Kc.^2-1);
a=(Kc/R).^2;

for i=1:length(w)
    f(i)=D0*alfa+D0*(1-alfa)*sum(a.*B*w(i)^2./(a.^2*D0^2+w(i)^2));
end


function g=BesselKernelsCylinder(n)
for l=1:n
    g(l)=fzero(inline('besselj(0,x)-besselj(1,x)/x'),2+(l-1)*pi);
end
function f=DwPlanar(w,d,D0,alfa,n)
% d = distance between planes
k=1:n;

B=8*d^2./(2*k-1).^4/pi^4;
a=(2*k-1).^2*pi^2/d^2;

for i=1:length(w)
    f(i)=D0*alfa+D0*(1-alfa)*sum(a.*B*w(i)^2./(a.^2*D0^2+w(i)^2));
end

f=f';
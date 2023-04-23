function[n]=LU_dec(an,bn,cn,fn)
% bn=diag(Matn);
% an=diag(Matn,-1);
% cn=diag(Matn,1);
Np=length(bn);
alphan(1) = bn(1);
for i=1:Np-1
    betan(i)=an(i)/alphan(i);
    alphan(i+1)=bn(i+1)-betan(i)*cn(i);
end
% Solution of Lv = f %    
vn(1) = fn(1);
for i = 2:Np
    vn(i) = fn(i) - betan(i-1)*vn(i-1);
end
% Solution of U*fi = v %    

tempn = vn(Np)/alphan(Np);
n(Np)=tempn;
for i = (Np-1):-1:1
    tempn = (vn(i)-cn(i)*n(i+1))/alphan(i);
    n(i) = tempn;
end
end
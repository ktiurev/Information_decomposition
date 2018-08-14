function [M1, M2, m1, m2, P_gen, str_gen] = Generation(C,x3,y3,n,N,Alphabet)

%GENERATION Summary of this function goes here% 

P_gen = zeros(n,n,n);
str_gen = '';
M1 = zeros(n,n);
M2 = zeros(n,n);
x = zeros(1,n);
x(x3) = 1.0;
y = zeros(1,n);
y(y3) = 1.0;


for i=1:N

% Generate P(Z) = sum_X(sum_Y(P(X,Y,Z)))    
        c2 = squeeze(C(x3,y3,:));
        z1 = c2/sum(c2);
        z2 = zeros(1,n); z2(1)=0; z2(1,2:n+1) = z1;
        z3 = sum(rand >= cumsum(z2));
        Z = Alphabet(z3);
        z=zeros(1,n); z(z3)=1;
        str_gen = [str_gen, ' ', Z, ' '];
        P_gen(x3,y3,z3) = P_gen(x3,y3,z3) + 1;
% Update matrices M1 and M2
        M1 = M1 + (z'*y);%(z.*y')
        M2 = M2 + (z'*x);%(z.*x')


% Generate P(X) = sum_Y(sum_Z(P(Y,Z,X)))    
        c2 = squeeze(C(y3,z3,:));
        x1 = c2/sum(c2);
        x2 = zeros(1,n); x2(1)=0; x2(1,2:n+1) = x1;
        x3 = sum(rand >= cumsum(x2));
        X = Alphabet(x3);
        x=zeros(1,n); x(x3)=1;
        str_gen = [str_gen, ' ', X,' '];
        P_gen(y3,z3,x3) = P_gen(y3,z3,x3) + 1;
% Update matrices M1 and M2
        M1 = M1 + (x'*z);%(x.*z');
        M2 = M2 + (x'*y);%(x.*y');


% Generate P(Y) = sum_X(sum_Z(P(Z,X,Y)))    
        c2 = squeeze(C(z3,x3,:));
        y1 = c2/sum(c2);
        y2 = zeros(1,n); y2(1)=0; y2(1,2:n+1) = y1;
        y3 = sum(rand >= cumsum(y2));
        Y = Alphabet(y3);
        y=zeros(1,n); y(y3)=1;
        str_gen = [str_gen, ' ', Y,' '];
        P_gen(z3,x3,y3) = P_gen(z3,x3,y3) + 1;
% Update matrices M1 and M2
        M1 = M1 + (y'.*x);%(y.*x');
        M2 = M2 + (y'.*z);%(y.*z');
end


% Column - normalize M1 and M2
m1 = M1;
m2 = M2;
pxy = M1/sum(sum(M1));
pxz = M2/sum(sum(M2));
M1 = pxy./sum(pxy,1);
M2 = pxz./sum(pxz,1);
P_gen = P_gen/sum(sum(sum(P_gen)));
    

end


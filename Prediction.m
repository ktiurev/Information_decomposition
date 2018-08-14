function [P_test, str_test] = Prediction(M1,M2,x3,y3,n,N,Alphabet)

P_test = zeros(n,n,n);
str_test = '';
x = zeros(1,n);
x(x3) = 1.0;
y = zeros(1,n);
y(y3) = 1.0;

    for i=1:N

        % Generate Z    
        c2 = M1*y' + M2*x';
        z1 = c2/sum(c2);
        z2 = zeros(1,n+1); z2(1)=0; z2(1,2:n+1) = z1;
        z3 = sum(rand >= cumsum(z2));
        Z = Alphabet(z3);
        z=zeros(1,n); z(z3)=1;
        str_test = [str_test, ' ', Z,' '];

        P_test(x3,y3,z3) = P_test(x3,y3,z3) + 1;

        % Generate X
        c2 = M1*z' + M2*y';
        x1 = c2/sum(c2);
        x2 = zeros(1,n+1); x2(1)=0; x2(1,2:n+1) = x1;
        x3 = sum(rand >= cumsum(x2));
        X = Alphabet(x3);
        x=zeros(1,n); x(x3)=1;
        str_test = [str_test, ' ', X,' '];

        P_test(y3,z3,x3) = P_test(y3,z3,x3) + 1;

        % Generate Y
        c2 = M1*x' + M2*z';
        y1 = c2/sum(c2);
        y2 = zeros(1,n+1); y2(1)=0; y2(1,2:n+1) = y1;
        y3 = sum(rand >= cumsum(y2));
        Y = Alphabet(y3);
        y=zeros(1,n); y(y3)=1;
        str_test = [str_test, ' ', Y,' '];

        P_test(z3,x3,y3) = P_test(z3,x3,y3) + 1;

    end
    
P_test = P_test/sum(sum(sum(P_test)));    
end


clear all

    % Number of generated symbols
    N = 10000;

    % Size and type of Alphabet
    n=4;
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    
% Build generative function

 % Construct joint probabilities
Pyz = rand(n,n);
Pyz = Pyz/sum(sum(Pyz));
Pxz = rand(n,n);
Pxz = Pxz/sum(sum(Pxz));

    % Construct magrinal probability of z
Pz = rand(n,1);
Pz = Pz/sum(Pz);

    % Construct conditional probabilities
Py_z = Pyz./Pz';
Px_z = Pxz./Pz';

    % Column-normalize conditional probabilities
Py_z = Py_z./sum(Py_z,1);
Px_z = Px_z./sum(Px_z,1);

    % Build full joint probability -- conditionally independent
for i=1:n
    for j=1:n
        for k=1:n
            Pxyz(i,j,k) = Px_z(i,k).*Py_z(j,k)*Pz(k);
        end
    end
end
    
    
    
    
 
        
% Fix probability distribution by hands
     Pxyz = rand(n,n,n);
     Pxyz(1,:,2) = 4.0;
     Pxyz(2,:,3) = 4.0;
     Pxyz(3,:,4) = 4.0;
     Pxyz(4,:,1) = 4.0;
     Pxyz = Pxyz/sum(sum(sum(Pxyz)));
    
    
% Check the probability distribution properties
% Check conditional independence
    Pxz_test = squeeze(sum(Pxyz,2));% --- Sum_y [P(x,y,z)]
    Pyz_test = squeeze(sum(Pxyz,1));% --- Sum_x [P(x,y,z)]
'Check the conditional independence'
    if( squeeze(max(max(Px_z.*Pz' - Pxz_test)))<10^-6 && squeeze(max(max(Py_z.*Pz' - Pyz_test)))<10^-6 ) 
        'ok'
    else
        'nope'
    end
'Check that the total probability is 1'
    if((sum(sum(sum(Pxyz)))-1.0)<10^-7) 
        'ok'
    else
        'nope'
    end

% Mutual information from M1 and M2

% Purely synergetic: XOR
%n=2;    
%p=zeros(n,n,n);
%p(1,1,1)=1;
%p(1,2,2)=1;
%p(2,1,2)=1;
%p(2,2,1)=1;
%p = p./sum(sum(sum(p)));
%Pxyz = p;


% Only unique from the first source
%clear all
%n=4;    
%p=rand(n,n,n);
%p(1,:,1) = 10;
%p(2,:,2) = 10;
%p(3,:,3) = 10;
%p(4,:,4) = 10;
%p=p./sum(sum(sum(p)));
%Pxyz = p;


% Highly redundant
%clear all
%n=4;    
%p=rand(n,n,n);
%p(1,1,1) = 10;
%p(2,2,2) = 10;
%p(3,3,3) = 10;
%p(4,4,4) = 10;
%p=p./sum(sum(sum(p)));
%Pxyz = p;




% Conditional independence: P(y,z|x) = P(y|x)P(z|x)
 if(1==0)
    Py_x = Py_x./sum(Py_x,1);
    Pz_x = Pz_x./sum(Pz_x,1);
	for i=1:n
        for j=1:n
            for k=1:n
                Pyzx(i,j,k) = Py_x(i,k)*Pz_x(j,k)*Px(k);
            end
        end
    end
	Pxyz = permute(Pyzx,[2 3 1]);
	Pxyz = Pxyz/sum(sum(sum(Pxyz)));
end


pxy = squeeze(sum(Pxyz,3));
pyx = permute(pxy,[2 1]);
pxz = squeeze(sum(Pxyz,2));
pzx = permute(pxz,[2 1]);
pyz = squeeze(sum(Pxyz,1));
pzy = permute(pyz,[2 1]);
px = squeeze(sum(sum(Pxyz,3),2));
py = squeeze(sum(sum(Pxyz,3),1));
pz = squeeze(sum(sum(Pxyz,2),1));


for i=1:n
    for j=1:n
        px_z(i,j) = pxz(i,j)/pz(j);
        py_z(i,j) = pyz(i,j)/pz(j);
        px_y(i,j) = pxy(i,j)/py(j);
        pz_y(i,j) = pzy(i,j)/py(j);
        py_x(i,j) = pyx(i,j)/px(j);
        pz_x(i,j) = pzx(i,j)/px(j);
    end
end

% Conditional mutual information I(x,y|Z)
pxyz = Pxyz;
pxzy = permute(Pxyz, [1 3 2]);
pyzx = permute(Pxyz, [2 3 1]);
for i=1:n
    for j=1:n
        for k=1:n
            pxy_z(i,j,k) = pxyz(i,j,k)/pz(k);
            pxz_y(i,j,k) = pxzy(i,j,k)/py(k);
            pyz_x(i,j,k) = pyzx(i,j,k)/px(k);
        end
    end
end

for i=1:n
    for j=1:n
        hxz(i,j) = pxz(i,j) .* log2(pxz(i,j) ./ (px(i) .* pz(j)));
        hxy(i,j) = pxy(i,j) .* log2(pxy(i,j) ./ (px(i) .* py(j)));
        hyz(i,j) = pyz(i,j) .* log2(pyz(i,j) ./ (py(i) .* pz(j)));
        for k=1:n
            hxy_z(i,j,k) = pxyz(i,j,k) .* log2(pxyz(i,j,k).*pz(k) ./ (pxz(i,k) .* pyz(j,k)));
            hxz_y(i,j,k) = pxzy(i,j,k) .* log2(pxzy(i,j,k).*py(k) ./ (pxy(i,k) .* pzy(j,k)));
            hyz_x(i,j,k) = pyzx(i,j,k) .* log2(pyzx(i,j,k).*px(k) ./ (pyx(i,k) .* pzx(j,k)));
        end
    end
    
end

hyz_x(isnan(hyz_x))=0.0;
hxz_y(isnan(hxz_y))=0.0;
hxy_z(isnan(hxy_z))=0.0;
hxz(isnan(hxz))=0.0;
hyz(isnan(hyz))=0.0;
hxy(isnan(hxy))=0.0;

 
'Conditional mutual information I(X,Y|Z) is'
Ixy_z = sum(sum(sum(hxy_z)))

'Conditional mutual information I(X,Z|Y) is'
Ixz_y = sum(sum(sum(hxz_y)))

'Conditional mutual information I(Y,Z|X) is'
Iyz_x = sum(sum(sum(hyz_x)))

'Mutual information I(X,Z) is'
Ixz = sum(sum(sum(hxz)))    

'Mutual information I(X,Y) is'
Ixy = sum(sum(sum(hxy)))

'Mutual information I(Y,Z) is'
Iyz = sum(sum(sum(hyz)))    
    
    

I_XYwithZ =  Iyz + Ixz_y
%I_XYwithZ = Ixz + Iyz_x



% Unique informaiton UI(Z,Y \ X)  
Pzy = squeeze(sum(Pxyz,1))'; Pzx = squeeze(sum(Pxyz,2))';
[UI2, Qp] = UI(Pzy, Pzx);
P = permute(Qp,[3 2 1]);
if( max(max(squeeze(sum(Pxyz,1))'-squeeze(sum(P,1))'))<10^-3 && max(max(squeeze(sum(Pxyz,2))'-squeeze(sum(P,2))'))<10^-3)
	'Correct margins'
end

% Unique informaiton UI(Z,X \ Y)
[UI1, Qp] = UI(Pzx, Pzy);
P = permute(Qp,[3 2 1]);
if( max(max(squeeze(sum(Pxyz,2))'-squeeze(sum(P,1))'))<10^-3 && max(max(squeeze(sum(Pxyz,1))'-squeeze(sum(P,2))'))<10^-3)
	'Correct margins'
end    
    
% Shared information 
Sh = Iyz - UI2
Sh = Ixz - UI1


% Synergetic information
Syn = I_XYwithZ - Sh - UI1 - UI2





% Generate first letters
    randvec = (1/n)*ones(1,n+1);
    randvec(1)=0;
    x0 = sum(rand >= cumsum(randvec));
    y0 = sum(rand >= cumsum(randvec));
    letter_x0 = Alphabet(x0);
    letter_y0 = Alphabet(y0);


% Generate the sequence with joint distribution P(X,Y,Z) starting from [x0, y0]
[M1, M2, m1, m2, P_gen, str_gen] = Generation(Pxyz,x0,y0,n,N,Alphabet);


% Generate the sequence with associative matrices M1 and M2 starting from [x0, y0]
[P_test, str_test] = Prediction(M1,M2,x0,y0,n,N,Alphabet);



   
% Sequence generated with joint probabilities P(X,Y,Z)
    textData = split(str_gen,newline);
    documents = tokenizedDocument(textData);
    bag = bagOfNgrams(documents,'NgramLengths',[2 3]);
    top = topkngrams(bag,10,'NGramLengths',3)

% Sequence generated with M1 and M2
    textData = split(str_test,newline);
    documents = tokenizedDocument(textData);
    bag2 = bagOfNgrams(documents,'NgramLengths',[2 3]);
    top2 = topkngrams(bag2,10,'NGramLengths',3)





    
    P_test = P_test + 10^-16;
    P_gen = P_gen + 10^-16;
    kl = P_gen.*log2(P_gen./P_test);
    sum(sum(sum(kl)));
    rep = 1;
    KulLeb(rep) = sum(sum(sum((kl))));
    'Kullback--Leibler divergence is'
    mean(KulLeb)
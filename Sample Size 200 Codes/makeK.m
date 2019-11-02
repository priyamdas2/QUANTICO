function P = makeK(M,difford)
if nargin<2
    difford=2;
end
P = eye(M);
if difford>0
    for d=1:difford
        P = diff(P);
    end
end
P = P'*P;
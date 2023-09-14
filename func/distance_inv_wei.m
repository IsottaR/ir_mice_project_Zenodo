function D=distance_inv_wei(W_)

n_=length(W_);
D=inf(n_);                                                  % distance matrix
D(1:n_+1:end)=0;

for u=1:n_
    S=true(1,n_);                                           % distance permanence (true is temporary)
    W1_=W_;
    V=u;
    while 1
        S(V)=0;                                             % distance u->V is now permanent
        W1_(:,V)=0;                                         % no in-edges as already shortest
        for v=V
            T=find(W1_(v,:));                               % neighbours of shortest nodes
            D(u,T)=min([D(u,T);D(u,v)+W1_(v,T)]);           % smallest of old/new path lengths
        end
        
        minD=min(D(u,S));
        if isempty(minD)||isinf(minD)                       % isempty: all nodes reached;
            break,                                          % isinf: some nodes cannot be reached
        end;
        
        V=find(D(u,:)==minD);
    end
end
D=1./D;                                                     % invert distance
D(1:n_+1:end)=0;
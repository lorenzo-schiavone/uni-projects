function [R, G, L, C, s] = componentMatrices_nonDC(types, vals, e)
% conductance
idx = find(types=='R');
v = 1./vals(idx);
G = sparse(idx,idx,v,e,e);
idx = find(types=='V' | types=='L');
v = ones(numel(idx),1);
G  = G + sparse(idx,idx, v, e, e);
% resistance
idx = find(types=='R');
v = ones(numel(idx),1);
R = sparse(idx,idx,v,e,e);
idx = find(types=='I' | types=='C');
v = ones(numel(idx),1);
R  = R + sparse(idx,idx, v, e, e);
% source
idx = find(types == 'V' | types == 'I');
v1 = vals(idx,1); v2 = vals(idx,2); v2(isnan(v2))=0;
v = v1.*exp(1i*v2);
s = sparse(idx,1,v,e,1);
% inductance
idx = find(types == 'L');
v = vals(idx);
L = sparse(idx, idx, v, e,e);
% capaticance
idx = find(types == 'C');
v = vals(idx);
C = sparse(idx,idx,v,e,e);
end
function [R, G, s] = componentMatricesDC(types, vals, e)
% conductance
idx = find(types=='R');
v = 1./vals(idx);
G = sparse(idx,idx,v,e,e);
idx = find(types=='V');
v = ones(numel(idx),1);
G  = G + sparse(idx,idx, v, e, e);
% resistance
idx = find(types=='R');
v = ones(numel(idx),1);
R = sparse(idx,idx,v,e,e);
idx = find(types=='I');
v = ones(numel(idx),1);
R  = R + sparse(idx,idx, v, e, e);
% source
idx = find(types == 'V' | types == 'I');
v = vals(idx);
s = sparse(idx,1,v,e,1);
end
function A = incidenceMatrix(nodes, n, e)
    cols = 1:e;
    Aplus = sparse(nodes(:,2), cols, ones(e,1), n, e);
    Aminus = sparse(nodes(:,1), cols, -ones(e,1), n, e);
    Ac = Aplus + Aminus;
    A = Ac(1:end-1,:);
end
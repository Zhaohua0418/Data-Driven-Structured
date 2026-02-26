function mi = info_transform(info)
    [m, T] = size(info.U);
    info.T = T;
    info.Xp = info.X(:,2:T+1);
    info.Xm = info.X(:,1:T);
    info.Um = info.U(:,1:T);
    mi = info;
end
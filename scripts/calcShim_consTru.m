function coef = calcShim_consTru(MapB0, MAP_Shims, currLimL, currLimU)

%% Apply consTru
coef = - MapB0*pinv(MAP_Shims) ;
S = svd(MAP_Shims);
idx = length(S);
while sum(coef>currLimU) || sum(coef<currLimL)
    coef = - MapB0*pinv(MAP_Shims,S(idx));
    idx = idx-1;
end
coef = transpose(coef);


function value = GetNddoRepulsionIntegral(obj, atomA, mu, nu, atomB, lambda, sigma)

if(atomA.atomType ~= 255 && atomB.atomType ~= 255)
    rAB = obj.molecule.GetDistanceAtoms(atomA, atomB);
elseif(atomA.atomType ~= 255 && atomB.atomType == 255)
    rAB = obj.molecule.GetDistanceAtomEpc(atomA, atomB);
elseif(atomA.atomType == 255 && atomB.atomType ~= 255)
    rAB = obj.molecule.GetDistanceAtomEpc(atomB, atomA);
else
    throw(MException('Mndo:GetNddoRepulsionIntegral', 'GetDistance method not implemented.'));
end

if(mu == 1 && nu == 1 && lambda == 1 && sigma == 1)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    
    % (29) in [DT_1977]
elseif(mu == 1 && nu == 1 && lambda == 4 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 2, rAB);
    value = temp1 + temp2;
    
elseif(mu == 1 && nu == 1 && lambda == 2 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 3, rAB);
    value = temp1 + temp2;
    
    % (30) in [DT_1977]
elseif(mu == 1 && nu == 1 && lambda == 3 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 4, rAB);
    value = temp1 + temp2;
    
    % (31) in [DT_1977]
elseif(mu == 4 && nu == 4 && lambda == 1 && sigma == 1)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 1, rAB);
    value = temp1 + temp2;
    
elseif(mu == 2 && nu == 2 && lambda == 1 && sigma == 1)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 1, rAB);
    value = temp1 + temp2;
    
    % (32) in [DT_1977]
elseif(mu == 3 && nu == 3 && lambda == 1 && sigma == 1)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 1, rAB);
    value = temp1 + temp2;
    
    % (33) in [DT_1977]
elseif(mu == 4 && nu == 4 && lambda == 4 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 2, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 2, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
elseif(mu == 2 && nu == 2 && lambda == 2 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 3, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 3, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
    % (34) in [DT_1977]
elseif(mu == 4 && nu == 4 && lambda == 2 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 3, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 3, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
elseif(mu == 2 && nu == 2 && lambda == 4 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 2, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 2, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
    % (35) in [DT_1977]
elseif(mu == 4 && nu == 4 && lambda == 3 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 4, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 4, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
elseif(mu == 2 && nu == 2 && lambda == 3 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 4, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 4, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
    % (36) in [DT_1977]
elseif(mu == 3 && nu == 3 && lambda == 4 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 2, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 2, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
elseif(mu == 3 && nu == 3 && lambda == 2 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 3, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 3, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
    % (37) in [DT_1977]
elseif(mu == 3 && nu == 3 && lambda == 3 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 4, rAB);
    temp3 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 1, rAB);
    temp4 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 4, rAB);
    value = temp1 + temp2 + temp3 + temp4;
    
    % (38) in [DT_1977]
elseif(mu == 1 && nu == 3 && lambda == 1 && sigma == 1)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 1, rAB);
    value = temp1;
    
elseif(mu == 3 && nu == 1 && lambda == 1 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
    % (39) in [DT_1977]
elseif(mu == 1 && nu == 3 && lambda == 4 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 2, rAB);
    value = temp1 + temp2;
    
elseif(mu == 3 && nu == 1 && lambda == 4 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 1 && nu == 3 && lambda == 2 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 3, rAB);
    value = temp1 + temp2;
    
elseif(mu == 3 && nu == 1 && lambda == 2 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
    % (40) in [DT_1977]
elseif(mu == 1 && nu == 3 && lambda == 3 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 1, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 4, rAB);
    value = temp1 + temp2;
    
elseif(mu == 3 && nu == 1 && lambda == 3 && sigma == 3)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
    % (41) in [DT_1977]
elseif(mu == 1 && nu == 1 && lambda == 1 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 10, rAB);
    value = temp1;
    
elseif(mu == 1 && nu == 1 && lambda == 3 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
    % (42) in [DT_1977]
elseif(mu == 4 && nu == 4 && lambda == 1 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 10, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 10, rAB);
    value = temp1 + temp2;
    
elseif(mu == 4 && nu == 4 && lambda == 3 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 2 && nu == 2 && lambda == 1 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 10, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 3, 10, rAB);
    value = temp1 + temp2;
    
elseif(mu == 2 && nu == 2 && lambda == 3 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
    % (43) in [DT_1977]
elseif(mu == 3 && nu == 3 && lambda == 1 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 1, 10, rAB);
    temp2 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 4, 10, rAB);
    value = temp1 + temp2;
    
elseif(mu == 3 && nu == 3 && lambda == 3 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
    % (44) in [DT_1977]
elseif(mu == 1 && nu == 4 && lambda == 1 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 8, 8, rAB);
    value = temp1;
    
elseif(mu == 4 && nu == 1 && lambda == 1 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 1 && nu == 4 && lambda == 4 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 4 && nu == 1 && lambda == 4 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
elseif(mu == 1 && nu == 2 && lambda == 1 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 9, 9, rAB);
    value = temp1;
    
elseif(mu == 2 && nu == 1 && lambda == 1 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 1 && nu == 2 && lambda == 2 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 2 && nu == 1 && lambda == 2 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
    % (45) in [DT_1977]
elseif(mu == 1 && nu == 3 && lambda == 1 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 10, rAB);
    value = temp1;
    
elseif(mu == 3 && nu == 1 && lambda == 1 && sigma == 3)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 1 && nu == 3 && lambda == 3 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 3 && nu == 1 && lambda == 3 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
    % (46) in [DT_1977]
elseif(mu == 1 && nu == 4 && lambda == 4 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 8, 5, rAB);
    value = temp1;
    
elseif(mu == 4 && nu == 1 && lambda == 4 && sigma == 3)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 1 && nu == 4 && lambda == 3 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 4 && nu == 1 && lambda == 3 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
elseif(mu == 1 && nu == 2 && lambda == 2 && sigma == 3)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 9, 6, rAB);
    value = temp1;
    
elseif(mu == 2 && nu == 1 && lambda == 2 && sigma == 3)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 1 && nu == 2 && lambda == 3 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 2 && nu == 1 && lambda == 3 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
    % (47) in [DT_1977]
elseif(mu == 4 && nu == 3 && lambda == 1 && sigma == 4)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 5, 8, rAB);
    value = temp1;
    
elseif(mu == 3 && nu == 4 && lambda == 1 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 4 && nu == 3 && lambda == 4 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 3 && nu == 4 && lambda == 4 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
elseif(mu == 2 && nu == 3 && lambda == 1 && sigma == 2)
    temp1 = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 6, 9, rAB);
    value = temp1;
    
elseif(mu == 3 && nu == 2 && lambda == 1 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 2 && nu == 3 && lambda == 2 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 3 && nu == 2 && lambda == 2 && sigma == 1)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
    % (48) in [DT_1977]
elseif(mu == 4 && nu == 3 && lambda == 4 && sigma == 3)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 5, 5, rAB);
    
elseif(mu == 3 && nu == 4 && lambda == 4 && sigma == 3)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 4 && nu == 3 && lambda == 3 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 3 && nu == 4 && lambda == 3 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
elseif(mu == 2 && nu == 3 && lambda == 2 && sigma == 3)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 6, 6, rAB);
    
elseif(mu == 3 && nu == 2 && lambda == 2 && sigma == 3)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 2 && nu == 3 && lambda == 3 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 3 && nu == 2 && lambda == 3 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
    % (49) in [DT_1977] and p19 in [MOPAC_1990]
elseif(mu == 4 && nu == 2 && lambda == 4 && sigma == 2)
    value = 0.5*(obj.GetNddoRepulsionIntegral(atomA, mu, mu, atomB, mu, mu)...
    -obj.GetNddoRepulsionIntegral(atomA, mu, mu, atomB, nu, nu));
    
elseif(mu == 2 && nu == 4 && lambda == 4 && sigma == 2)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
    
elseif(mu == 4 && nu == 2 && lambda == 2 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
    
elseif(mu == 2 && nu == 4 && lambda == 2 && sigma == 4)
    value = obj.GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
    
    % d-orbitals
elseif(mu == 5 || mu == 6 || mu == 7 || mu == 8 || mu == 9 ||...
        nu == 5 || nu == 6 || nu == 7 || nu == 8 || nu == 9 ||...
        lambda == 5 || lambda == 6 || lambda == 7 || lambda  == 8 || lambda == 9 ||...
        sigma == 5 || sigma == 6 || sigma == 7 || sigma  == 8 || sigma == 9)
    
    throw(MException('Mndo:GetNddoRepulsionIntegral', 'd-orbitals not implemented.'));
    
else
    value = 0.0;
end
end

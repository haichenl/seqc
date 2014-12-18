function value = Cpp_GetNddoRepulsionIntegral(obj, atomA, mu, nu, atomB, lambda, sigma)
% if(atomA.atomType ~= 255 && atomB.atomType ~= 255)
    rAB = obj.molecule.GetDistanceAtoms(atomA, atomB);
% elseif(atomA.atomType ~= 255 && atomB.atomType == 255)
%     rAB = obj.molecule.GetDistanceAtomEpc(atomA, atomB);
% elseif(atomA.atomType == 255 && atomB.atomType ~= 255)
%     rAB = obj.molecule.GetDistanceAtomEpc(atomB, atomA);
% else
%     throw(MException('Mndo:GetNddoRepulsionIntegral', 'GetDistance method not implemented.'));
% end
DAvec = obj.AtomGetNddoDerivedParameterDVec(atomA);
DBvec = obj.AtomGetNddoDerivedParameterDVec(atomB);
% if(atomA.atomType ~= 255 && atomB.atomType ~= 255)
    rhoAvec = obj.AtomGetNddoDerivedParameterRhoVec(atomA);
    rhoBvec = obj.AtomGetNddoDerivedParameterRhoVec(atomB);
% elseif(atomA.atomType ~= 255 && atomB.atomType == 255)
%     rhoAvec = zeros(3, 1);
%     rhoBvec = zeros(3, 1);
% elseif(atomA.atomType == 255 && atomB.atomType ~= 255)
%     rhoAvec = zeros(3, 1);
%     rhoBvec = zeros(3, 1);
% else
%     throw(MException('Mndo:GetSemiEmpiricalMultipoleInteractionl', 'GetNddoDerivedParameterRho method not implemented.'));
% end

value = SEQC.Mndo.Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, lambda, sigma, rAB);

end
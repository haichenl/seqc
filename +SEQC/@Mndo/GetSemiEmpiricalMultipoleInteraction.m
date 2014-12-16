function value = GetSemiEmpiricalMultipoleInteraction(obj, atomA, atomB, multipoleA, multipoleB, rAB)
DA = obj.AtomGetNddoDerivedParameterD(atomA, multipoleA);
DB = obj.AtomGetNddoDerivedParameterD(atomB, multipoleB);
if(atomA.atomType ~= 255 && atomB.atomType ~= 255)
    rhoA = obj.AtomGetNddoDerivedParameterRho(atomA, multipoleA);
    rhoB = obj.AtomGetNddoDerivedParameterRho(atomB, multipoleB);
elseif(atomA.atomType ~= 255 && atomB.atomType == 255)
    rhoA = 0.0;
    rhoB = 0.0;
elseif(atomA.atomType == 255 && atomB.atomType ~= 255)
    rhoA = 0.0;
    rhoB = 0.0;
else
    throw(MException('Mndo:GetSemiEmpiricalMultipoleInteractionl', 'GetNddoDerivedParameterRho method not implemented.'));
end
a = rhoA + rhoB;

% Eq. (52) in [DT_1977]
if(multipoleA == 1 && multipoleB == 1)
    value = 1.0/sqrt(rAB*rAB + a*a);
    
    % Eq. (53) in [DT_1977]
elseif(multipoleA == 1 && multipoleB == 10)
    temp1 = ((rAB+DB)*(rAB+DB)) + (a*a);
    temp2 = ((rAB-DB)*(rAB-DB)) + (a*a);
    value = 1.0/sqrt(temp1)/2.0 - 1.0/sqrt(temp2)/2.0;
    
elseif(multipoleA == 10 && multipoleB == 1)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    value = value * (-1.0);
    
    % Eq. (54) in [DT_1977]
elseif(multipoleA == 1 && multipoleB == 2)
    temp1 = (rAB*rAB) + (4.0*DB*DB) + (a*a);
    temp2 = (rAB*rAB) + (a*a);
    value = 1.0/sqrt(temp1)/2.0 - 1.0/sqrt(temp2)/2.0;
    
elseif(multipoleA == 2 && multipoleB == 1)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    
elseif(multipoleA == 1 && multipoleB == 3)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, multipoleA, 2, rAB);
    
elseif(multipoleA == 3 && multipoleB == 1)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    
    % Eq. (55) in [DT_1977]
elseif(multipoleA == 1 && multipoleB == 4)
    temp1 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (a*a);
    temp2 = (rAB*rAB) + (a*a);
    temp3 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (a*a);
    value = 1.0/sqrt(temp1)/4.0 - 1.0/sqrt(temp2)/2.0 + 1.0/sqrt(temp3)/4.0;
    
elseif(multipoleA == 4 && multipoleB == 1)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    
    % Eq. (56) in [DT_1977]
elseif(multipoleA == 8 && multipoleB == 8)
    temp1 = (rAB*rAB) + ((DA-DB)*(DA-DB)) + (a*a);
    temp2 = (rAB*rAB) + ((DA+DB)*(DA+DB)) + (a*a);
    value = 1.0/sqrt(temp1)/2.0 - 1.0/sqrt(temp2)/2.0;
    
elseif(multipoleA == 9 && multipoleB == 9)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 8, 8, rAB);
    
    % Eq. (57) in [DT_1977]
elseif(multipoleA == 10 && multipoleB == 10)
    temp1 = ((rAB+DA-DB)*(rAB+DA-DB)) + (a*a);
    temp2 = ((rAB+DA+DB)*(rAB+DA+DB)) + (a*a);
    temp3 = ((rAB-DA-DB)*(rAB-DA-DB)) + (a*a);
    temp4 = ((rAB-DA+DB)*(rAB-DA+DB)) + (a*a);
    value = 1.0/sqrt(temp1)/4.0 - 1.0/sqrt(temp2)/4.0...
        -1.0/sqrt(temp3)/4.0 + 1.0/sqrt(temp4)/4.0;
    
    % Eq. (58) in [DT_1977]
elseif(multipoleA == 8 && multipoleB == 5)
    temp1 = ((rAB-DB)*(rAB-DB)) + ((DA-DB)*(DA-DB)) + (a*a);
    temp2 = ((rAB-DB)*(rAB-DB)) + ((DA+DB)*(DA+DB)) + (a*a);
    temp3 = ((rAB+DB)*(rAB+DB)) + ((DA-DB)*(DA-DB)) + (a*a);
    temp4 = ((rAB+DB)*(rAB+DB)) + ((DA+DB)*(DA+DB)) + (a*a);
    value =-1.0/sqrt(temp1)/4.0 + 1.0/sqrt(temp2)/4.0...
        +1.0/sqrt(temp3)/4.0 - 1.0/sqrt(temp4)/4.0;
    
elseif(multipoleA == 5 && multipoleB == 8)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    value = value * (-1.0);
    
elseif(multipoleA == 9 && multipoleB == 6)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 8, 5, rAB);
    
elseif(multipoleA == 6 && multipoleB == 9)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    value = value * (-1.0);
    
    % Eq. (59) in [DT_1977]
elseif(multipoleA == 10 && multipoleB == 2)
    temp1 = ((rAB+DA)*(rAB+DA)) + (4.0*DB*DB) + (a*a);
    temp2 = ((rAB-DA)*(rAB-DA)) + (4.0*DB*DB) + (a*a);
    temp3 = ((rAB+DA)*(rAB+DA)) + (a*a);
    temp4 = ((rAB-DA)*(rAB-DA)) + (a*a);
    value =-1.0/sqrt(temp1)/4.0 + 1.0/sqrt(temp2)/4.0...
        +1.0/sqrt(temp3)/4.0 - 1.0/sqrt(temp4)/4.0;
    
elseif(multipoleA == 2 && multipoleB == 10)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    value = value * (-1.0);
    
elseif(multipoleA == 10 && multipoleB == 3)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 10, 2, rAB);
    
elseif(multipoleA == 3 && multipoleB == 10)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    value = value * (-1.0);
    
    % Eq. (60) in [DT_1977]
elseif(multipoleA == 10 && multipoleB == 4)
    temp1 = ((rAB+DA-2.0*DB)*(rAB+DA-2.0*DB)) + (a*a);
    temp2 = ((rAB-DA-2.0*DB)*(rAB-DA-2.0*DB)) + (a*a);
    temp3 = ((rAB+DA+2.0*DB)*(rAB+DA+2.0*DB)) + (a*a);
    temp4 = ((rAB-DA+2.0*DB)*(rAB-DA+2.0*DB)) + (a*a);
    temp5 = ((rAB+DA)*(rAB+DA)) + (a*a);
    temp6 = ((rAB-DA)*(rAB-DA)) + (a*a);
    value =-1.0/sqrt(temp1)/8.0 + 1.0/sqrt(temp2)/8.0...
        -1.0/sqrt(temp3)/8.0 + 1.0/sqrt(temp4)/8.0...
        +1.0/sqrt(temp5)/4.0 - 1.0/sqrt(temp6)/4.0;
    
elseif(multipoleA == 4 && multipoleB == 10)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    value = value * (-1.0);
    
    % Eq. (61) in [DT_1977]
elseif(multipoleA == 2 && multipoleB == 2)
    temp1 = (rAB*rAB) + 4.0*((DA-DB)*(DA-DB)) + (a*a);
    temp2 = (rAB*rAB) + 4.0*((DA+DB)*(DA+DB)) + (a*a);
    temp3 = (rAB*rAB) + (4.0*DA*DA) + (a*a);
    temp4 = (rAB*rAB) + (4.0*DB*DB) + (a*a);
    temp5 = (rAB*rAB) + (a*a);
    value = 1.0/sqrt(temp1)/8.0 + 1.0/sqrt(temp2)/8.0...
        -1.0/sqrt(temp3)/4.0 - 1.0/sqrt(temp4)/4.0...
        +1.0/sqrt(temp5)/4.0;
    
elseif(multipoleA == 3 && multipoleB == 3)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, 2, rAB);
    
    % Eq. (62) in [DT_1977]
elseif(multipoleA == 2 && multipoleB == 3)
    temp1 = (rAB*rAB) + (4.0*DA*DA) + (4.0*DB*DB)+ (a*a);
    temp2 = (rAB*rAB) + (4.0*DA*DA) + (a*a);
    temp3 = (rAB*rAB) + (4.0*DB*DB) + (a*a);
    temp4 = (rAB*rAB) + (a*a);
    value = 1.0/sqrt(temp1)/4.0 - 1.0/sqrt(temp2)/4.0...
        -1.0/sqrt(temp3)/4.0 + 1.0/sqrt(temp4)/4.0;
    
elseif(multipoleA == 3 && multipoleB == 2)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    
    % Eq. (63) in [DT_1977]
elseif(multipoleA == 2 && multipoleB == 4)
    temp1 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (4.0*DA*DA) + (a*a);
    temp2 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (4.0*DA*DA) + (a*a);
    temp3 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (a*a);
    temp4 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (a*a);
    temp5 = (rAB*rAB) + (4.0*DA*DA) + (a*a);
    temp6 = (rAB*rAB) + (a*a);
    value = 1.0/sqrt(temp1)/8.0 + 1.0/sqrt(temp2)/8.0...
        -1.0/sqrt(temp3)/8.0 - 1.0/sqrt(temp4)/8.0...
        -1.0/sqrt(temp5)/4.0 + 1.0/sqrt(temp6)/4.0;
    
elseif(multipoleA == 4 && multipoleB == 2)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    
elseif(multipoleA == 3 && multipoleB == 4)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 2, multipoleB, rAB);
    
elseif(multipoleA == 4 && multipoleB == 3)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
    
    % Eq. (64) in [DT_1977]
elseif(multipoleA == 4 && multipoleB == 4)
    temp1 = ((rAB+2.0*DA-2.0*DB)*(rAB+2.0*DA-2.0*DB)) + (a*a);
    temp2 = ((rAB+2.0*DA+2.0*DB)*(rAB+2.0*DA+2.0*DB)) + (a*a);
    temp3 = ((rAB-2.0*DA-2.0*DB)*(rAB-2.0*DA-2.0*DB)) + (a*a);
    temp4 = ((rAB-2.0*DA+2.0*DB)*(rAB-2.0*DA+2.0*DB)) + (a*a);
    temp5 = ((rAB+2.0*DA)*(rAB+2.0*DA)) + (a*a);
    temp6 = ((rAB-2.0*DA)*(rAB-2.0*DA)) + (a*a);
    temp7 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (a*a);
    temp8 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (a*a);
    temp9 = (rAB*rAB) + (a*a);
    value = 1.0/sqrt(temp1)/16.0 + 1.0/sqrt(temp2)/16.0...
        +1.0/sqrt(temp3)/16.0 + 1.0/sqrt(temp4)/16.0...
        -1.0/sqrt(temp5)/8.0 - 1.0/sqrt(temp6)/8.0...
        -1.0/sqrt(temp7)/8.0 - 1.0/sqrt(temp8)/8.0...
        +1.0/sqrt(temp9)/4.0;
    
    % Eq. (65) in [DT_1977]
elseif(multipoleA == 5 && multipoleB == 5)
    temp1 = ((rAB+DA-DB)*(rAB+DA-DB)) + ((DA-DB)*(DA-DB)) + (a*a);
    temp2 = ((rAB+DA-DB)*(rAB+DA-DB)) + ((DA+DB)*(DA+DB)) + (a*a);
    temp3 = ((rAB+DA+DB)*(rAB+DA+DB)) + ((DA-DB)*(DA-DB)) + (a*a);
    temp4 = ((rAB+DA+DB)*(rAB+DA+DB)) + ((DA+DB)*(DA+DB)) + (a*a);
    temp5 = ((rAB-DA-DB)*(rAB-DA-DB)) + ((DA-DB)*(DA-DB)) + (a*a);
    temp6 = ((rAB-DA-DB)*(rAB-DA-DB)) + ((DA+DB)*(DA+DB)) + (a*a);
    temp7 = ((rAB-DA+DB)*(rAB-DA+DB)) + ((DA-DB)*(DA-DB)) + (a*a);
    temp8 = ((rAB-DA+DB)*(rAB-DA+DB)) + ((DA+DB)*(DA+DB)) + (a*a);
    value = 1.0/sqrt(temp1)/8.0 - 1.0/sqrt(temp2)/8.0...
        -1.0/sqrt(temp3)/8.0 + 1.0/sqrt(temp4)/8.0...
        -1.0/sqrt(temp5)/8.0 + 1.0/sqrt(temp6)/8.0...
        +1.0/sqrt(temp7)/8.0 - 1.0/sqrt(temp8)/8.0;
    
elseif(multipoleA == 6 && multipoleB == 6)
    value = obj.GetSemiEmpiricalMultipoleInteraction(atomA, atomB, 5, 5, rAB);
    
    % Eq. (66) in [DT_1977]
elseif(multipoleA == 7 && multipoleB == 7)
    temp1 = (rAB*rAB) + 2.0*((DA-DB)*(DA-DB)) + (a*a);
    temp2 = (rAB*rAB) + 2.0*((DA+DB)*(DA+DB)) + (a*a);
    temp3 = (rAB*rAB) + 2.0*(DA*DA) + 2.0*(DB*DB) + (a*a);
    value = 1.0/sqrt(temp1)/4.0 + 1.0/sqrt(temp2)/4.0...
        -1.0/sqrt(temp3)/2.0;
    
else
    throw(MException('Mndo:GetSemiEmpiricalMultipoleInteractionl', 'Multipole type wrong.'));
end
end

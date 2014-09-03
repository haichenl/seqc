classdef (Abstract) Atom < handle
    
    properties (SetAccess = protected)
        
        xyz;
        atomType;
        atomicMass;
        valence;
        realSphericalHarmonicsIndices;
        valenceShellType;
        firstAOIndex;
        numberValenceElectrons;
        vdWCoefficient;
        vdWRadii;
        
        coreCharge;
        
        % singleton holding parameters
        paramPool;
        
    end
    
    properties (SetAccess = private)
        index;
    end
    
    methods
        
        function obj = Atom(ind)
            obj.index = ind;
        end
        
        function SetXyz(obj, xyz_in)
            if(prod(size(xyz_in)==[3,1]))
                % do nothing
            elseif(prod(size(xyz_in)==[1,3]))
                xyz_in = xyz_in';
            else
                throw(MException('Atom:SetXyz', 'Input xyz vector dimension wrong.'));
            end
            obj.xyz = xyz_in;
        end
        
        function res = GetCoreMass(obj)
            res = obj.atomicMass - obj.numberValenceElectrons;
        end
        
        function res = GetValenceSize(obj)
            res = length(obj.valence);
        end
        
        function SetCoreCharge(obj, charge)
            obj.coreCharge = charge;
        end
        
        function SetFirstAOIndex(obj, firstAOIndex_in)
            obj.firstAOIndex = firstAOIndex_in;
        end
        
        function res = GetLastAOIndex(obj)
            res = obj.firstAOIndex + length(obj.valence) - 1;
        end
        
        function res = GetEffectivePrincipalQuantumNumber(~, shellType) % ShellType shellType
%             if(shellType == 1) % EnumShell.kShell)
%                 res = 1.0;
%             elseif(shellType == 2) % EnumShell.lShell)
%                 res = 2.0;
%             elseif(shellType == 3) % EnumShell.mShell)
%                 res = 3.0;
%             elseif(shellType == 4) % EnumShell.nShell)
%                 res = 3.7;
%             else
%                 throw(MException('Atom:GetEffectivePrincipalQuantumNumber', 'Shell type wrong.'));
%             end
            res = double(shellType);
            for i = 1:length(res)
                if(res(i) == 1 || res(i) == 2 || res(i) == 3)
                    % do nothing
                elseif(res(i) == 4)
                    res(i) = 3.7;
                else
                    throw(MException('Atom:GetEffectivePrincipalQuantumNumber', 'Shell type wrong.'));
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function res = GetRealAngularPartAO(~, theta, phi, orbital) % OrbitalType orbital
            switch(orbital)
                case 1 % EnumOrbital.s
                    res = power(4.0*pi,-0.5);
                case 2 % EnumOrbital.py
                    res = power(3.0/(4.0*pi),0.5)*sin(theta)*sin(phi);
                case 3 % EnumOrbital.pz
                    res = power(3.0/(4.0*pi),0.5)*cos(theta);
                case 4 % EnumOrbital.px
                    res = power(3.0/(4.0*pi),0.5)*sin(theta)*cos(phi);
                case 5 % EnumOrbital.dxy
                    res = power(15.0/(16.0*pi),0.5)*power(sin(theta),2.0)*sin(2.0*phi);
                case 6 % EnumOrbital.dyz
                    res = power(15.0/(16.0*pi),0.5)*sin(2.0*theta)*sin(phi);
                case 7 % EnumOrbital.dzz
                    res = power(5.0/(16.0*pi),0.5)*(3.0*power(cos(theta),2.0) - 1.0);
                case 8 % EnumOrbital.dzx
                    res = power(15.0/(16.0*pi),0.5)*sin(2.0*theta)*cos(phi);
                case 9 % EnumOrbital.dxxyy
                    res = power(15.0/(16.0*pi),0.5)*power(sin(theta),2.0)*cos(2.0*phi);
                otherwise
                    throw(MException('Atom:GetRealAngularPartAO', 'Orbital type not possible'));
            end
        end
        
        function res = GetRadialPartAO(~, dr, orbitalExponent, shell) % ShellType shell
            principalQuantumNumber = double(shell);
            temp1 = power(2.0*orbitalExponent,principalQuantumNumber+0.5);
            temp2 = power(factorial(2*principalQuantumNumber),-0.5);
            res = temp1*temp2*power(dr,principalQuantumNumber-1)*exp(-1.0*orbitalExponent*dr);
        end
        
    end
    
end

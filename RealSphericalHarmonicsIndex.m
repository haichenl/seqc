classdef RealSphericalHarmonicsIndex < handle
    
    properties (SetAccess = private)
        l;
        m;
    end
    
    methods
        
        function obj = RealSphericalHarmonicsIndex(arg1, arg2) % OrbitalType orbitalType
            if(nargin == 1)
                orbitalType = arg1;
                if(orbitalType == OrbitalType.s)
                    obj.l = 0;
                    obj.m = 0;
                elseif(orbitalType == OrbitalType.py)
                    obj.l = 1;
                    obj.m = -1;
                elseif(orbitalType == OrbitalType.pz)
                    obj.l = 1;
                    obj.m = 0;
                elseif(orbitalType == OrbitalType.px)
                    obj.l = 1;
                    obj.m = 1;
                elseif(orbitalType == OrbitalType.dxy)
                    obj.l = 2;
                    obj.m = -2;
                elseif(orbitalType == OrbitalType.dyz)
                    obj.l = 2;
                    obj.m = -1;
                elseif(orbitalType == OrbitalType.dzz)
                    obj.l = 2;
                    obj.m = 0;
                elseif(orbitalType == OrbitalType.dzx)
                    obj.l = 2;
                    obj.m = 1;
                elseif(orbitalType == OrbitalType.dxxyy)
                    obj.l = 2;
                    obj.m = 2;
                else
                    throw(MException('RealSphericalHarmonicsIndex:RealSphericalHarmonicsIndex', 'Orbital type wrong.'));
                end
            elseif(nargin == 2)
                obj.l = arg1;
                obj.m = arg2;
            else
                throw(MException('RealSphericalHarmonicsIndex:RealSphericalHarmonicsIndex', 'Input argument wrong'));
            end
        end
        
        function res = GetL(obj)
            res = obj.l;
        end
        
        function res = GetM(obj)
            res = obj.m;
        end
        
    end
    
end
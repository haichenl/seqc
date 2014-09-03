classdef RealSphericalHarmonicsIndex < handle
    
    properties (SetAccess = private)
        l;
        m;
    end
    
    methods
        
        function obj = RealSphericalHarmonicsIndex(arg1, arg2) % OrbitalType orbitalType
            if(nargin == 1)
                orbitalType = arg1;
                if(orbitalType == 1)
                    obj.l = 0;
                    obj.m = 0;
                elseif(orbitalType == 2)
                    obj.l = 1;
                    obj.m = -1;
                elseif(orbitalType == 3)
                    obj.l = 1;
                    obj.m = 0;
                elseif(orbitalType == 4)
                    obj.l = 1;
                    obj.m = 1;
                elseif(orbitalType == 5)
                    obj.l = 2;
                    obj.m = -2;
                elseif(orbitalType == 6)
                    obj.l = 2;
                    obj.m = -1;
                elseif(orbitalType == 7)
                    obj.l = 2;
                    obj.m = 0;
                elseif(orbitalType == 8)
                    obj.l = 2;
                    obj.m = 1;
                elseif(orbitalType == 9)
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
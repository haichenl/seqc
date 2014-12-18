classdef (Sealed) ParamPoolManagerPm3 < FITSEQC.ParamPoolManagerBase
    
    methods (Access = protected)
        
        function res = ParamsLengthOfAPool(~, paramPool)
            res = length(paramPool.pm3ValidParams);
        end
        
        function res = GetIthPoolParams(obj, i)
            res = obj.paramPoolsVec{i}.GetPm3Params();
        end
        
        function SetIthPoolParams(obj, i, ithPoolParams)
            obj.paramPoolsVec{i}.SetPm3Params(ithPoolParams);
        end
        
        % private constructor
        function obj = ParamPoolManagerPm3()
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = FITSEQC.ParamPoolManagerPm3();
            end
            singleObj = localObj;
        end
        
    end
    
end
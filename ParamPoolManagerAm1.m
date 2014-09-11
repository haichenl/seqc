classdef (Sealed) ParamPoolManagerAm1 < ParamPoolManagerBase
    
    methods (Access = protected)
        
        function res = ParamsLengthOfAPool(~, paramPool)
            res = length(paramPool.am1ValidParams);
        end
        
        function res = GetIthPoolParams(obj, i)
            res = obj.paramPoolsVec{i}.GetAm1Params();
        end
        
        function SetIthPoolParams(obj, i, ithPoolParams)
            obj.paramPoolsVec{i}.SetAm1Params(ithPoolParams);
        end
        
        % private constructor
        function obj = ParamPoolManagerAm1()
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = ParamPoolManagerAm1();
            end
            singleObj = localObj;
        end
        
    end
    
end
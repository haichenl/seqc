classdef (Sealed) ParamPoolManagerCndo2 < FITSEQC.ParamPoolManagerBase
    
    methods (Access = protected)
        
        function res = ParamsLengthOfAPool(~, paramPool)
            res = length(paramPool.cndo2ValidParams);
        end
        
        function res = GetIthPoolParams(obj, i)
            res = obj.paramPoolsVec{i}.GetCndo2Params();
        end
        
        function SetIthPoolParams(obj, i, ithPoolParams)
            obj.paramPoolsVec{i}.SetCndo2Params(ithPoolParams);
        end
        
        % private constructor
        function obj = ParamPoolManagerCndo2()
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = FITSEQC.ParamPoolManagerCndo2();
            end
            singleObj = localObj;
        end
        
    end
    
end
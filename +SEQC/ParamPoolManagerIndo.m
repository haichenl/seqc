classdef (Sealed) ParamPoolManagerIndo < SEQC.ParamPoolManagerBase
    
    methods (Access = protected)
        
        function res = ParamsLengthOfAPool(~, paramPool)
            res = length(paramPool.indoValidParams);
        end
        
        function res = GetIthPoolParams(obj, i)
            res = obj.paramPoolsVec{i}.GetIndoParams();
        end
        
        function SetIthPoolParams(obj, i, ithPoolParams)
            obj.paramPoolsVec{i}.SetIndoParams(ithPoolParams);
        end
        
        % private constructor
        function obj = ParamPoolManagerIndo()
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = SEQC.ParamPoolManagerIndo();
            end
            singleObj = localObj;
        end
        
    end
    
end
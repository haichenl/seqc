classdef (Abstract) ParamPoolManagerBase < handle
    
    properties (SetAccess = private)
        
        paramPoolsVec = {};
        paramStartIndices = [];
        paramEndIndices = [];
        
    end
    
    methods (Access = public)
        
        function AddParamPool(obj, paramPool)
            paramsLengthOfAPool = obj.ParamsLengthOfAPool(paramPool);
            if(isempty(obj.paramPoolsVec))
                obj.paramPoolsVec{1} = paramPool;
                obj.paramStartIndices(1) = 1;
                obj.paramEndIndices(1) = obj.paramStartIndices(1) + paramsLengthOfAPool - 1;
                return;
            end
            for i = 1:length(obj.paramPoolsVec)
                if(obj.paramPoolsVec{i} == paramPool)
                    return;
                end
            end
            obj.paramPoolsVec{end+1} = paramPool;
            obj.paramStartIndices(end+1) = obj.paramEndIndices(end) + 1;
            obj.paramEndIndices(end+1) = obj.paramStartIndices(end) + paramsLengthOfAPool - 1;
        end
        
        function res = GetAllParams(obj)
            res = zeros(obj.paramEndIndices(end), 1);
            for i = 1:length(obj.paramPoolsVec)
                res(obj.paramStartIndices(i):obj.paramEndIndices(i)) = obj.GetIthPoolParams(i);
            end
        end
        
        function SetAllParams(obj, allParams)
            for i = 1:length(obj.paramPoolsVec)
                obj.paramPoolsVec{i}.SetCndo2Params(allParams(obj.paramStartIndices(i):obj.paramEndIndices(i)));
            end
        end
        
    end
    
    methods (Abstract, Access = protected)
        
        res = ParamsLengthOfAPool(obj, paramPool)
        
        res = GetIthPoolParams(obj, i)
        
        SetIthPoolParams(obj, i, ithPoolParams)

    end
    
end
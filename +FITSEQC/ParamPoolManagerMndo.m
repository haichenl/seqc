classdef (Sealed) ParamPoolManagerMndo < FITSEQC.ParamPoolManagerBase
    
    methods (Access = protected)
        
        function res = ParamsLengthOfAPool(~, paramPool)
            res = length(paramPool.mndoValidParams);
        end
        
        function res = GetIthPoolParams(obj, i)
            res = obj.paramPoolsVec{i}.GetMndoParams();
        end
        
        function SetIthPoolParams(obj, i, ithPoolParams)
            obj.paramPoolsVec{i}.SetMndoParams(ithPoolParams);
        end
        
        % private constructor
        function obj = ParamPoolManagerMndo()
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = FITSEQC.ParamPoolManagerMndo();
            end
            singleObj = localObj;
        end
        
    end
    
end
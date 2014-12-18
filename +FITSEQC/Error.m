classdef Error < handle
    
    properties
        
        weightEnergyTotal = 1;
        weightEnergyHOMO = 1;
        weightEnergyLUMO = 1;
        weightEnergyEnvironment = 1;
        
    end
    
    properties (Access = private)
        
        paramManager;
        
        lowLevelCell;
        highLevelCell;
        
        referenceVector;
        indicesNotIncludingReferences;
        
        lowLevelEnergyPool;
        highLevelEnergyPool;
        
    end
    
    properties (Access = private, Constant)
        
        TotalNumberOfTypes = 4;
        IndexEnergyTotal = 1;
        IndexEnergyHOMO = 2;
        IndexEnergyLUMO = 3;
        IndexEnergyEnvironment = 4;
        
    end
    
    methods
        
        function obj = Error(input)
            lowLevelCell_ = input.lowLevelCell;
            highLevelCell_ = input.highLevelCell;
            refNumbers = input.refNumbers;
            paramManager_ = input.paramManager;
            
            if(length(lowLevelCell_) ~= length(highLevelCell_))
                throw(MException('Error:Error','low/high level numbers do not agree.'));
            end
            obj.lowLevelCell = lowLevelCell_;
            obj.highLevelCell = highLevelCell_;
            
            obj.lowLevelEnergyPool = zeros(length(obj.lowLevelCell), 4);
            obj.highLevelEnergyPool = zeros(length(obj.highLevelCell), 4);
            for iHigh = 1:length(obj.highLevelCell)
                obj.highLevelEnergyPool(iHigh, :) = ...
                    obj.EnergyColumnFromOneModel(obj.highLevelCell{iHigh});
            end
            
            obj.referenceVector = zeros(1, length(obj.lowLevelCell));
            refNumbers = [reshape(refNumbers,1,[]) length(obj.lowLevelCell)+1];
            obj.indicesNotIncludingReferences = [];
            for iref = 1:length(refNumbers)-1
                obj.referenceVector(refNumbers(iref):refNumbers(iref+1)-1) ...
                    = refNumbers(iref).*ones(1, refNumbers(iref+1)-refNumbers(iref));
                tempIndices = refNumbers(iref)+1:refNumbers(iref+1)-1;
                obj.indicesNotIncludingReferences(end+1:end+length(tempIndices)) ...
                    = tempIndices;
            end
            
            obj.paramManager = paramManager_;
        end
        
        function currParams = CurrentParameters(obj)
            currParams = obj.paramManager.GetAllParams();
        end
        
        function errorColumn = ErrorFunction(obj, newParameters)
            obj.paramManager.SetAllParams(newParameters);
            obj.UpdateLowLevelEnergyPool();
            errorColumn = [ ...
                obj.RelativeErrorColumn(obj.IndexEnergyTotal) .* obj.weightEnergyTotal; ...
                obj.AbsoluteErrorColumn(obj.IndexEnergyHOMO) .* obj.weightEnergyHOMO;
                obj.AbsoluteErrorColumn(obj.IndexEnergyLUMO) .* obj.weightEnergyLUMO;
                obj.AbsoluteErrorColumn(obj.IndexEnergyEnvironment) .* obj.weightEnergyEnvironment;
                ];
        end
        
    end
    
    methods (Access = private)
        
        function column = EnergyColumnFromOneModel(obj, model)
            column = zeros(obj.TotalNumberOfTypes, 1);
            column(obj.IndexEnergyTotal) = model.energyTotal;
            column(obj.IndexEnergyHOMO) = model.energyHOMO;
            column(obj.IndexEnergyLUMO) = model.energyLUMO;
            column(obj.IndexEnergyEnvironment) = model.energyEnvironment;
        end
        
        function UpdateLowLevelEnergyPool(obj)
            for iLow = 1:length(obj.lowLevelCell)
                obj.lowLevelCell{iLow}.DoSCF();
                obj.lowLevelEnergyPool(iLow, :) = ...
                    obj.EnergyColumnFromOneModel(obj.lowLevelCell{iLow});
            end
        end
        
        function column = AbsoluteErrorColumn(obj, index)
            column = obj.lowLevelEnergyPool(:, index) ...
                - obj.highLevelEnergyPool(:, index);
        end
        
        function column = RelativeErrorColumn(obj, index)
            columnLow = obj.lowLevelEnergyPool(:, index);
            columnHigh = obj.highLevelEnergyPool(:, index);
            column = (columnLow - columnHigh) ...
                - (columnLow(obj.referenceVector) - columnHigh(obj.referenceVector));
            column = column(obj.indicesNotIncludingReferences);
        end
        
    end
    
end
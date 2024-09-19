

classdef PropConfig
    properties
        EnergyInt = false;           % No additional structure options
        PhaseInt = false;            % Integration of the phase
        
        EnrgyProj = false;           % Projection on the osculating Hamiltonian
        OrthoProj = false;           % Orthogonal projection (preservation of the bilinear function)

        PhaseCnt = false;            % Use of the phase for kinematic control

        SundmanTransformation = "1"; % Sundman transformation  
    end

    methods
        function [obj] = PropConfig()
            obj.EnergyInt = false;
        end
    end
end
function mainFunction()
    runSim();

    function runSim()
    
        model = createpde(1);

        % Define the block geometry
        width = 256e-9;   % Width 
        len = 256e-9;     % Length
        height = 22e-9;   % Height 
        block = multicuboid(width, len, height);
        model.Geometry = block;

        % Generate the mesh for the geometry
        generateMesh(model, 'Hmax', min([width, len, height])/10);

        % Define material properties and coefficients
        epsilon_r_dielectric = 7.5; % Silicon Nitride
        specifyCoefficients(model, 'm', 0, 'd', 0, 'c', @(location,state)materialProperty(location, epsilon_r_dielectric, width, len, height), 'a', 0, 'f', 0);

        % Apply boundary conditions
       
        applyBoundaryCondition(model, 'dirichlet', 'Face', 1:model.Geometry.NumFaces, 'u', 0); %default to entire block gnd

        % Solve the PDE
        result = solvepde(model);

        % Plot the solution
        figure;
        pdeplot3D(model, 'ColorMapData', result.NodalSolution);
        title('Electric Potential');
    end

    function c = materialProperty(location, epsilon_dielectric, width, len, height)
        % Initialize permittivity with dielectric value
        c = epsilon_dielectric * ones(size(location.x));

        % Mimic grounding and applying voltage by adjusting permittivity
        % bottom wire region -up to 2e-9 meters
        groundRegionZmax = 2e-9;

        % Voltage-applied region: top down to height - 2e-9 meters
        voltageRegionZmin = height - 2e-9;

       
        highEpsilon = 1e6; % high permit. for conductor

        % Identify locations within the specified regions
        inGroundRegion = location.z <= groundRegionZmax;
        inVoltageRegion = location.z >= voltageRegionZmin;

        % Assign conductor permit. to regions
        c(inGroundRegion) = highEpsilon;
        c(inVoltageRegion) = highEpsilon;
    end
end

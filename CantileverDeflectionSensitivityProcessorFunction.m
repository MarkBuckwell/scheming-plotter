function [ DeflArray , SensLocation ] = CantileverDeflectionSensitivityProcessorFunction ( )

%% Imports and processes sets of .txt AFM force-distance data to calculate the deflection sensitivity.
% Data cropped to middle third between minimum and maximum voltage (y-axis).
% This assumes the region will be a straight line, corresponding to the
% contact regime. The line is fitted using a first order polynomial and the
% mean, standard deviation and standard error extracted.

%% Selects and imports data.
addpath ( cd ) ;
[ FlipArray , DeflSensVals , DeflkVals , DeflDataHeaders , NFC , aH , aV , rH , rV ] = CantileverDeflectionFileSelector ( ) ;

%% Assumes extend deflection (VEx) and height sensor (nmEx),
% are the first and penultimate data columns, respectively. Retract colums
% are second and final for deflection and height sensor, respectively.

iFirstRun = 0 ;

while iFirstRun < 1

    FirstFit = questdlg ( 'Trial fitting attempts, final or close program?' , 'Trial fitting attempts, final or close program?' , 'Trials' , 'Final' , 'Close' , 'Trials' ) ;

    switch FirstFit  % User option to maximise figure window for ease of inspection.

        case 'Trials'

            TrialFits = 1 ;

        case 'Final'

            TrialFits = 0 ;

        case 'Close' 

            TrialFits = 1 ;

    end

    %% Check for sampling parameters file
    FileName = strcat ( 'samplingparameters.txt' ) ;
    if exist ( FileName , 'file' ) > 0 % Checks to see if sampling parameters file exists.

        ParametersFile = fopen ( FileName ) ; % If present, parameters are opened.
        FileParameters = textscan ( ParametersFile , '%f' ) ;
        Parameters = FileParameters { 1 } ;
        CropStartEx = Parameters ( 1 ) ; % Crop start point for extend data.
        CropSizeEx = Parameters ( 2 ) ; % Crop size for extend data.
        CropStartRt = Parameters ( 3 ) ; % Crop start point for retract data.
        CropSizeRt = Parameters ( 4 ) ; % Crop size for retract data.

    else % Otherwise, use basic/standard set of parameters.

        %     CropSize = ( round ( ( Samples / 15 ) , 0 ) - 10 ) ; % Number of
        % samples for cropped data, if dependent on total number of samples.
        CropStartEx = 7 ; % Crop start point for extend data.
        CropSizeEx = 35 ; % Crop size for extend data.
        CropStartRt = 7 ; % Crop start point for retract data.
        CropSizeRt = 35 ; % Crop size for retract data.

     end

    %% With final sampling parameters enterd, the data can be fit.
    if TrialFits == 0

        CropArrayEx = zeros ( CropSizeEx + 1 , 2 , NFC ) ; % Array to contain cropped extend data.
        CropArrayRt = zeros ( CropSizeRt + 1 , 2 , NFC ) ; % Array to contain cropped retract data.
        SensArray = zeros ( NFC , 2 ) ; % Array to contain fitted deflection sensitivity values.

        figure ( 1 ) ;
        clf

        for iNFC = 1 : NFC

            CropArrayEx ( : , 1 , iNFC ) = FlipArray ( ( end - CropSizeEx - CropStartEx ) : ( end - CropStartEx ) , aH , iNFC ) ; % Crops cantilever extend height data.
            CropArrayRt ( : , 1 , iNFC ) = FlipArray ( ( end - CropSizeRt - CropStartRt ) : ( end - CropStartRt ) , rH , iNFC ) ; % Crops cantilever retract height data.
            CropArrayEx ( : , 2 , iNFC ) = FlipArray ( ( end - CropSizeEx - CropStartEx ) : ( end - CropStartEx ) , aV , iNFC ) ; % Crops cantilever extend voltage data.
            CropArrayRt ( : , 2 , iNFC ) = FlipArray ( ( end - CropSizeRt - CropStartRt ) : ( end - CropStartRt ) , rV , iNFC ) ; % Crops cantilever retract voltage data.

            LinearFitEx = polyfit ( CropArrayEx ( : , 1 , iNFC ) , CropArrayEx ( : , 2 , iNFC ) , 1 ) ; % Linear fit to extend data.
            LinearFitRt = polyfit ( CropArrayRt ( : , 1 , iNFC ) , CropArrayRt ( : , 2 , iNFC ) , 1 ) ; % Linear fit to retract data.

            SensArray ( iNFC , 1 ) = 1 / LinearFitEx ( 1 ) ;
            SensArray ( iNFC , 2 ) = 1 / LinearFitRt ( 1 ) ;

            X1 = CropArrayEx ( : , 1 , iNFC ) - min ( CropArrayEx ( : , 1 , iNFC ) ) ;
            Y1 = CropArrayEx ( : , 2 , iNFC ) - min ( CropArrayEx ( : , 2 , iNFC ) )  ;
            X2 = CropArrayRt ( : , 1 , iNFC ) - min ( CropArrayRt ( : , 1 , iNFC ) ) ;
            Y2 = CropArrayRt ( : , 2 , iNFC ) - min ( CropArrayRt ( : , 2 , iNFC ) ) ;

            Y3 = ( LinearFitEx ( 1 ) .* X1 ) ;
            Y4 = ( LinearFitRt ( 1 ) .* X2 ) ;

            hold on
            subplot ( 2 , 1 , 1 ) ;
            plot (  X1 , Y1 , 'o' ,  X2 , Y2 , 'o' , X1 , Y3 , X2 , Y4 ) ;
            Leg1 = strcat ( 'File' , { ' ' } , num2str ( iNFC ) ) ;
            Leg2 = strcat ( 'CropSizeEx' , { ' ' } , num2str ( CropSizeEx ) ) ;
            Leg3 = strcat ( 'CropStartEx' , { ' ' } , num2str ( CropStartEx ) ) ;
            Leg4 = strcat ( 'CropSizeRt' , { ' ' } , num2str ( CropSizeRt ) ) ;
            Leg5 = strcat ( 'CropStartRt' , { ' ' } , num2str ( CropStartRt ) ) ;
            title ( Leg1 { 1 } ) ;
            legend ( Leg2 { 1 } , Leg3 { 1 } , Leg4 { 1 } , Leg5 { 1 } , 'Location' , 'SouthEast' ) ;
            %waitforbuttonpress ;

        end

        savefig ( 'deflplot.fig' ) ; % Saves current fit, overwriting previous.


        %% Creates array for values and easy copying out of Matlab.

        DeflArray = cell ( 9 , 2 ) ;

        MeanSensEx = mean ( SensArray ( : , 1 ) ) ;
        MeanSensRt = mean ( SensArray ( : , 2 ) ) ;
        StDevSensEx = std ( SensArray ( : , 1 ) ) ;
        StDevSensRt = std ( SensArray ( : , 2 ) ) ;
        StErrorSensEx = StDevSensEx / sqrt ( NFC ) ;
        StErrorSensRt = StDevSensRt / sqrt ( NFC ) ;

        DeflArray { 1 } = 'MeanSensEx (nm/V)' ;
        DeflArray { 2 } = 'MeanSensRt (nm/V)' ;
        DeflArray { 3 } = 'MeanSens (nm/V)' ;
        DeflArray { 4 } = 'StDevSensEx (nm/V)' ;
        DeflArray { 5 } = 'StDevSensRt (nm/V)' ;
        DeflArray { 6 } = 'StDevSens (nm/V)' ;
        DeflArray { 7 } = 'StErrorSensEx (nm/V)' ;
        DeflArray { 8 } = 'StErrorSensRt (nm/V)' ;
        DeflArray { 9 } = 'StErrorSens (nm/V)' ;
        DeflArray { 10 } = MeanSensEx ;
        DeflArray { 11 } = MeanSensRt ;
        DeflArray { 12 } = ( MeanSensEx + MeanSensRt ) / 2 ;
        DeflArray { 13 } = StDevSensEx ;
        DeflArray { 14 } = StDevSensRt ;
        DeflArray { 15 } = ( StDevSensEx + StDevSensRt ) / 2 ;
        DeflArray { 16 } = StErrorSensEx ;
        DeflArray { 17 } = StErrorSensRt ;
        DeflArray { 18 } = ( StErrorSensEx + StErrorSensRt ) / 2 ;
        
        for i = 1 : NFC % Plotted to overlay all data as a final check that no spectra look anomalous.

            hold on

            subplot ( 2 , 1 , 2 ) ;
            plot ( FlipArray ( : , 1 , i ) , FlipArray ( : , 3 , i ) , 'o' ) ;
            plot ( FlipArray ( : , 2 , i ) , FlipArray ( : , 4 , i ) ) ;

        end
        
        ModulusArray = zeros ( NFC , 2 ) ;
        
        for i = 1 : NFC % Loop to calculate differences between approach and retract for each sample for energy dissipation data.
           
            [ Min , MinIndex ] = min ( FlipArray ( : , 3 , i ) ) ;
            ModulusArray ( i , 1 ) = sum ( FlipArray ( MinIndex : end , 3 , i ) - FlipArray ( MinIndex : end , 4 , i ) ) * ...
            ( max ( FlipArray ( : , 1 , i ) ) - FlipArray ( MinIndex , 1 , i ) ) ;
            ModulusArray ( i , 2 ) = max ( max ( vertcat ( FlipArray ( : , 3 , i ) , FlipArray ( : , 4 , i ) ) ) ) ;
            
        end
        
        iFirstRun = 1 ;

    %% Varied parameters to check optimum sampling of data.
    elseif TrialFits == 1

        CropSizeMat = ( 5 : 1 : 100 ) ; % Number of samples for cropped data.
        CropStartMat = ( 0 : 1 : 50 ) ; % Number of samples from maximum deflection value to begin sampling data for crop.
        ExArrayAll = zeros ( length ( CropSizeMat ) , length ( CropStartMat ) , NFC ) ; % Array to contain all fitted extend deflection sensitivity values.
        RtArrayAll = zeros ( length ( CropSizeMat ) , length ( CropStartMat ) , NFC ) ; % Array to contain all fitted retract deflection sensitivity values.
        ExArrayMean = zeros ( length ( CropSizeMat ) , length ( CropStartMat )  ) ; % Array to contain mean fitted extend deflection sensitivity values.
        RtArrayMean = zeros ( length ( CropSizeMat ) , length ( CropStartMat )  ) ; % Array to contain mean fitted retract deflection sensitivity values.
        ExArrayStDev = zeros ( length ( CropSizeMat ) , length ( CropStartMat )  ) ; % Array to contain fitted extend deflection sensitivity values standard deviations.
        RtArrayStDev = zeros ( length ( CropSizeMat ) , length ( CropStartMat )  ) ; % Array to contain fitted retract deflection sensitivity values standard deviations.

        for iStart = 1 : length ( CropStartMat )

            CropStart = CropStartMat ( iStart ) ;

            for iSize = 1 : length ( CropSizeMat )

                CropSize = CropSizeMat ( iSize ) ;
                CropArray = zeros ( CropSize + 1 , 4 , NFC ) ; % Array to contain cropped data.

                for iNFC = 1 : NFC

                    CropArray ( : , 1 , iNFC ) = FlipArray ( ( end - CropSize - CropStart ) : ( end - CropStart ) , 1 , iNFC ) ; % Cuts cantilever approach from height data, up to minimum voltage value.
                    CropArray ( : , 2 , iNFC ) = FlipArray ( ( end - CropSize - CropStart ) : ( end - CropStart ) , 2 , iNFC ) ; % Cuts cantilever approach from voltage data, up to minimum voltage value.
                    CropArray ( : , 3 , iNFC ) = FlipArray ( ( end - CropSize - CropStart ) : ( end - CropStart ) , 3 , iNFC ) ; % Cuts cantilever retract from height data, after minimum voltage value.
                    CropArray ( : , 4 , iNFC ) = FlipArray ( ( end - CropSize - CropStart ) : ( end - CropStart ) , 4 , iNFC ) ; % Cuts cantilever retract from voltage data, after minimum voltage value.

                    LinearFitEx = polyfit ( CropArray ( : , 1 , iNFC ) , CropArray ( : , 3 , iNFC ) , 1 ) ; % Linear fit to extend data.
                    LinearFitRt = polyfit ( CropArray ( : , 2 , iNFC ) , CropArray ( : , 4 , iNFC ) , 1 ) ; % Linear fit to retract data.

                    ExArrayAll ( iSize , iStart , iNFC ) = 1 / LinearFitEx ( 1 ) ;
                    RtArrayAll ( iSize , iStart , iNFC ) = 1 / LinearFitRt ( 1 ) ;

                end

                ExArrayMean ( iSize , iStart ) = mean ( ExArrayAll ( iSize , iStart , : ) ) ;
                ExArrayStDev ( iSize , iStart ) = std ( ExArrayAll ( iSize , iStart , : ) ) ;
                RtArrayMean ( iSize , iStart ) = mean ( RtArrayAll ( iSize , iStart , : ) ) ;
                RtArrayStDev ( iSize , iStart ) = std ( RtArrayAll ( iSize , iStart , : ) ) ;

            end

        end

        DiffArray = abs ( ExArrayMean - RtArrayMean ) ; % Absolute difference between extend and retract mean data.
        [ MinSizeAll , iSizeAll ] = min ( DiffArray ) ; % Minimum difference values and corresponding indices.
        [ MinSize , mStart  ] = min ( MinSizeAll ) ; % Minimum difference value and corresponding minimum difference CropStart values.
        mSize = iSizeAll ( mStart ) ; % Index of minimum difference CropSize value.

        MeanRange = 3 ; % To crop data scale for easier viewing.
    %     ExRange = [ ( mean ( mean ( ExArray ) ) - ( 2 * MeanRange ) ) , ( mean ( mean ( ExArray ) ) + MeanRange ) ] ;
    %     RtRange = [ ( mean ( mean ( RtArray ) ) - ( 2 * MeanRange ) ) , ( mean ( mean ( RtArray ) ) + MeanRange ) ] ;

        ExRange = [ ( min ( min ( ExArrayMean ) ) ) , ( min ( min ( ExArrayMean ) ) + MeanRange ) ] ;
        RtRange = [ ( min ( min ( RtArrayMean ) ) ) , ( min ( min ( RtArrayMean ) ) + MeanRange ) ] ;


        %% Plot mean data.

        figure ( 1 ) ;
        clf

        subplot ( 2 , 1 , 1 ) ;
        pcolor ( CropStartMat , CropSizeMat , ExArrayMean ) ;
        title ( 'Mean extend sensitivity (nm/V)' ) ;
        xlabel ( 'CropStart' ) ;
        ylabel ( 'CropSize' ) ;
        shading flat ;
        colorbar ;
%         caxis ( ExRange ) ;

        subplot ( 2 , 1 , 2 ) ;
        pcolor ( CropStartMat , CropSizeMat , RtArrayMean ) ;
        title ( 'Mean retract sensitivity (nm/V)' ) ;
        xlabel ( 'CropStart' ) ;
        ylabel ( 'CropSize' ) ;
        shading flat ;
        colorbar ;
%         caxis ( RtRange ) ;

        colormap ( jet ) ;
        set ( gcf , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

        %% Plot standard deviation data.

        figure ( 2 ) ;
        clf

        subplot ( 2 , 2 , 1 ) ;
        pcolor ( CropStartMat , CropSizeMat , ExArrayStDev ) ;
        title ( 'Extend sensitivity standard deviation (nm/V)' ) ;
        xlabel ( 'CropStart' ) ;
        ylabel ( 'CropSize' ) ;
        shading flat ;
        colorbar ;

        subplot ( 2 , 2 , 2 ) ;
        pcolor ( CropStartMat , CropSizeMat , RtArrayStDev ) ;
        title ( 'Retract sensitivity standard deviation (nm/V)' ) ;
        xlabel ( 'CropStart' ) ;
        ylabel ( 'CropSize' ) ;
        shading flat ;
        colorbar ;

        subplot ( 2 , 2 , [ 3 , 4 ] ) ;
        pcolor ( CropStartMat , CropSizeMat , DiffArray ) ;
        title ( 'Absolute mean sensitivity difference (nm/V)' ) ;
        xlabel ( strcat ( 'CropStart; min value' , { ' ' }  , num2str ( CropStartMat (  mStart ) ) ) ) ;
        ylabel ( strcat ( 'CropSize; min value' , { ' ' } , num2str ( CropSizeMat ( mSize ) ) ) ) ;
        shading flat ;
        colorbar ;

        colormap ( jet ) ;
        set ( gcf , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

        %% User input to choose sampling parameters, given test parameters in figures.

        figure ( 1 ) ;
        iSwitch = 0 ;
        iData = 0 ;

        while iSwitch < 2

            Accept = MFquestdlg ( [ 0.7 0.73 ] , 'Enter parameters or change figure?' , 'Enter parameters or change figure?' , 'Enter parameters' , ...
                    'Switch figure' , 'Change data range' , 'Enter parameters' ) ;

            switch Accept

            case 'Enter parameters' % Allow user to enter sampling parameters.
            % User can close dialog box to exit program here.

                UserPrompt = { 'CropStart extend value' , 'CropSize extend value' , 'CropStart retract value' , 'CropSize retract value' } ;
                PromptTitle = 'Enter sampling parameters' ;
                PromptLines = 1 ;
                DefaultAnswers = { ( num2str ( CropStartEx ) ) ( num2str ( CropSizeEx ) ) ( num2str ( CropStartRt ) ) ( num2str ( CropSizeRt ) ) } ;
                AlterProperties = inputdlg ( UserPrompt , PromptTitle , PromptLines , DefaultAnswers ) ;

                CropStartExtendUser = str2double ( AlterProperties { 1 } ) ; % Crop start point for extend data.
                CropSizeExtendUser = str2double ( AlterProperties { 2 } ) ; % Crop size for extend data.
                CropStartRetractUser = str2double ( AlterProperties { 3 } ) ; % Crop start point for retract data.
                CropSizeRetractUser = str2double ( AlterProperties { 4 } ) ; % Crop size for retract data.

                save ( 'samplingparameters.txt' , 'CropStartExtendUser', 'CropSizeExtendUser' , 'CropStartRetractUser' , 'CropSizeRetractUser' , '-ascii' ) ;

                iSwitch = 3 ;

            case 'Switch figure'

                if iSwitch == 0

                    figure ( 2 ) ;
                    iSwitch = 1 ;

                elseif iSwitch == 1

                    figure ( 1 ) ;
                    iSwitch = 0 ;

                end

            case 'Change data range'

                UserPrompt = { 'Upper range extend' , 'Lower range extend' , 'Upper range retract' , 'Lower range retract' } ;
                PromptTitle = 'Enter data range' ;
                PromptLines = 1 ;

                if iData == 0

                    DefaultAnswers = { ( num2str ( MeanRange ) ) '0' ( num2str ( MeanRange ) ) '0' } ;

                elseif iData == 1

                    DefaultAnswers = { ( num2str ( DataUpperEx ) ) ( num2str ( DataLowerEx ) ) ( num2str ( DataUpperRt ) ) ( num2str ( DataLowerRt ) ) } ;

                end

                AlterProperties = inputdlg ( UserPrompt , PromptTitle , PromptLines , DefaultAnswers ) ;  

                DataUpperEx = str2double ( AlterProperties { 1 } ) ;
                DataLowerEx = str2double ( AlterProperties { 2 } ) ;
                DataUpperRt = str2double ( AlterProperties { 3 } ) ;
                DataLowerRt = str2double ( AlterProperties { 4 } ) ;

                ExRange = [ ( min ( min ( ExArrayMean ) ) - DataLowerEx ) , ( min ( min ( ExArrayMean ) ) + DataUpperEx ) ] ;
                RtRange = [ ( min ( min ( RtArrayMean ) ) - DataLowerRt ) , ( min ( min ( RtArrayMean ) ) + DataUpperRt ) ] ;

                %% Replot mean data.

                figure ( 1 ) ;
                clf

                subplot ( 2 , 1 , 1 ) ;
                pcolor ( CropStartMat , CropSizeMat , ExArrayMean ) ;
                title ( 'Mean extend sensitivity (nm/V)' ) ;
                xlabel ( 'CropStart' ) ;
                ylabel ( 'CropSize' ) ;
                shading flat ;
                colorbar ;
                caxis ( ExRange ) ;

                subplot ( 2 , 1 , 2 ) ;
                pcolor ( CropStartMat , CropSizeMat , RtArrayMean ) ;
                title ( 'Mean retract sensitivity (nm/V)' ) ;
                xlabel ( 'CropStart' ) ;
                ylabel ( 'CropSize' ) ;
                shading flat ;
                colorbar ;
                caxis ( RtRange ) ;

                colormap ( jet ) ;
                set ( gcf , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

                %% Plot standard deviation data.

                figure ( 2 ) ;
                clf

                subplot ( 2 , 2 , 1 ) ;
                pcolor ( CropStartMat , CropSizeMat , ExArrayStDev ) ;
                title ( 'Extend sensitivity standard deviation (nm/V)' ) ;
                xlabel ( 'CropStart' ) ;
                ylabel ( 'CropSize' ) ;
                shading flat ;
                colorbar ;

                subplot ( 2 , 2 , 2 ) ;
                pcolor ( CropStartMat , CropSizeMat , RtArrayStDev ) ;
                title ( 'Retract sensitivity standard deviation (nm/V)' ) ;
                xlabel ( 'CropStart' ) ;
                ylabel ( 'CropSize' ) ;
                shading flat ;
                colorbar ;

                subplot ( 2 , 2 , [ 3 , 4 ] ) ;
                pcolor ( CropStartMat , CropSizeMat , DiffArray ) ;
                title ( 'Absolute mean sensitivity difference (nm/V)' ) ;
                xlabel ( strcat ( 'CropStart; min value' , { ' ' }  , num2str ( CropStartMat (  mStart ) ) ) ) ;
                ylabel ( strcat ( 'CropSize; min value' , { ' ' } , num2str ( CropSizeMat ( mSize ) ) ) ) ;
                shading flat ;
                colorbar ;

                colormap ( jet ) ;
                set ( gcf , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

                iData = 1 ;
                figure ( 1 ) ;

            end

        end

    end
    
end

SensLocation = cd ;

end

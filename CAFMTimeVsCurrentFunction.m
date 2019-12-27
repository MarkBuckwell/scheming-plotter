function [ OutputData , CharArray , SetThreshold , VIterationRange , h ] = CAFMTimeVsCurrentFunction ( )

set ( 0 , 'DefaultFigureWindowStyle' , 'Normal' ) ;

%% Imports and plots current data in time from  2-column AFM data.
% Exported as a .txt file from a 2D image file, should be a 0D scan in
% which the probe has been fixed at a point and the current sampled in time
% rather than space. Thus, trace and retrace current will be interpolated
% line by line into a full current measurement. All files should have the
% same number of headerlines and the same data columns, arranged in the
% same order.

addpath ( cd ) ;

[ FileGroup , DataPath ] = uigetfile ( '*.txt' , 'DialogTitle' ,...
    'Select files:' , 'MultiSelect' , 'on' ) ; % Gets file names and location.

% Used to account for the case where only a single file is selected.
SingleFile = double ( ischar ( FileGroup ) ) ;

if  SingleFile > 0

    NFC = 1 ;
    FileSet = strcat ( DataPath , FileGroup ) ;

else

    NFC = length ( FileGroup ) ;  % Number of files to import.
    FileSet = repmat ( { '' } , 1 , NFC ) ; % Generates cell array to place filenames into.

%% Generates array of files for analysis.

    for i = 1 : NFC

        % Concatenates path and file strings and adds to output array.
        FileSet ( i ) = strcat ( DataPath , FileGroup ( i ) ) ;

    end

end

cd ( DataPath ) ; % Sets target directory to current folder.

%% Recipient arrays.

CurrentAndTime = zeros ( 1 , 2 ) ; % Recipient matrix for current and time data.
VoltageAndTime = zeros ( 1 , 2 ) ; % Recipient matrix for resistance and time data.
SampleTime = zeros ( NFC , 3 )  ; % Recipient matrix for time data.
DataArray = repmat ( { '' } , 1 , 8 , NFC ) ; % Sets up recipient array for all data.
MissingSamples = repmat ( { ' ' } , NFC , size ( DataArray , 2 ) ) ; 
CurrentPerFile = repmat ( { [ 0 , 0 ] } , 1 , NFC ) ;
VoltagePerFile = CurrentPerFile ;

%% Import data.

for i = 1 : NFC

    if NFC == 1 

        FileName = char ( FileSet ) ; % Choose single file.
        disp ( FileGroup ) ;

    else

        FileName = char ( FileSet ( : , i ) ) ; % Choose file from set.
        disp ( FileGroup { i } ) ;

    end

    %% Get scan calibration values.
    FileChoice = fopen ( FileName , 'r' ) ;
    CalArray = textscan ( FileChoice, '%s' , 'Delimiter' , '\n'  ) ;
    CalArray = CalArray { 1 } ;

    % Header row and number of columns.
    HeaderRow = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
    'C-AFM_Current(V)' ) ) ) ) * 1 ) ; % Specifies the row containing headers.
    DataHeaders = textscan ( CalArray { HeaderRow } , '%s' , 'Delimiter' , ')' ) ;
    DataHeaders = DataHeaders { 1 } .' ;
    NumData = numel ( DataHeaders ) ;

    % Specifies import format; NumData cells containing data columns.
    FormatSpec = repmat ( '%f' , 1 , NumData ) ;

    % Height sensitivity.
    HeightSensCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'Sens. ZsensSens' ) ) ) ) ) ;
    HeightSens = textscan ( CalArray { HeightSensCell ( 1 ) } , '%s' ) ;
    HeightSens = str2double ( cell2mat ( ( HeightSens { 1 } ( 4 ) ) ) ) ;

    % Current sensitivity.
    CurrentSensCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'AFM Current Sensitivity' ) ) ) ) ) ;
    CurrentSens = textscan ( CalArray { CurrentSensCell ( 1 ) } , '%s' ) ;
    CurrentSens = str2double ( cell2mat ( ( CurrentSens { 1 } ( 4 ) ) ) ) ;

    % Applied voltage (in case raw voltage data has been full fitted).
    ApplVCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'DC Sample Bias' ) ) ) ) ) ;
    ApplV = textscan ( CalArray { ApplVCell ( 1 ) } , '%s' ) ;
    ApplV = abs ( str2double ( cell2mat ( ( ApplV { 1 } ( 9 ) ) ) ) ) ;

    % Samples per line.
    SampleFreqCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'Samps/line' ) ) ) ) ) ;
    SamplesPerLine = textscan ( CalArray { SampleFreqCell ( 1 ) } , '%s' ) ;
    SamplesPerLine = cell2mat ( ( SamplesPerLine { 1 } ( 2 ) ) ) ;
    SamplesPerLine ( end ) = [] ;
    SamplesPerLine = str2double ( SamplesPerLine ) ;

    % Number of scan lines.
    NumLinesCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'Number of lines' ) ) ) ) ) ;
    NumLines = textscan ( CalArray { NumLinesCell ( 2 ) } , '%s' ) ;
    NumLines = cell2mat ( ( NumLines { 1 } ( 4 ) ) ) ;
    NumLines ( end ) = [] ;
    NumLines = str2double ( NumLines ) ;

    % Scan rate (i.e. rate to perform both a trace and a retrace).
    RateCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'Scan Rate' ) ) ) ) ) ;
    ScanRate = textscan ( CalArray { RateCell ( 2 ) } , '%s' ) ;
    ScanRate = cell2mat ( ( ScanRate { 1 } ( 3 ) ) ) ;
    ScanRate ( end ) = [] ;
    ScanRate = str2double ( ScanRate ) ;

    % Time per sample.
    SampleTime ( i , 1 ) = 1 ./ ( ScanRate .* 2 .* SamplesPerLine ) ;

    % Opens and reads selected file data to temporary array.
    frewind ( FileChoice ) ;
    DataArray ( : , : , i ) = textscan ( FileChoice , FormatSpec , 'Delimiter' ,...
        '\t' , 'EmptyValue' , NaN , 'HeaderLines' , HeaderRow  , 'ReturnOnError' , false ) ;
    fclose ( FileChoice ) ;

    % Check that the data is the correct length.
    for j = 1 : size ( DataArray , 2 )

        RealLines ( i , j ) = size ( DataArray { 1 , j , i } , 1 ) ;

        if RealLines ( i , j ) ~= NumLines * SamplesPerLine

            disp ( strcat ( 'Data length is incorrect. There should be' ,...
                { ' ' } , num2str ( NumLines ) , ', but there are only' ,...
                { ' ' } , num2str ( RealLines ( i , j ) / SamplesPerLine ) ,...
                '. Data has been automatically corrected.' ) ) ;

            MissingSamples { i , j }  = zeros ( ( NumLines * SamplesPerLine ) -...
            RealLines ( i , j ), 1 ) + DataArray { 1 , j , i } ( end ) ;
            DataArray { 1 , j , i } = vertcat ( DataArray { 1 , j , i } ,...
                MissingSamples { i , j } ) ;

        end

    end

    %% Get indices of current and voltage data columns.
    iC = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
        'C-AFM' ) ) ) ) ) ;
    iV = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
        'Input' ) ) ) ) ) ;

    % Makes current absolute.
    DataArray { : , iC ( 1 ) , i } = abs ( DataArray { : , iC ( 1 ) , i } ) .* CurrentSens ;

    if size ( iC , 2 ) == 2

        DataArray { : , iC ( 2 ) , i } = abs ( DataArray { : , iC ( 2 ) , i } ) .* CurrentSens ;

    end

    %% Section to check data direction.
    % Scan direction - data must be flipped for up scans, as Nanoscpe
    % exports under the assumption that a 2D image has been acquired, not a
    % 0D image wherein the first line of a new image should correspond to
    % the top rather than the bottom. i.e. we are interested in the
    % evolution of current/voltage in time, not in space.

    DirCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'direction' ) ) ) ) ) ;
    ScanDir = textscan ( CalArray { DirCell ( 2 ) } , '%s' ) ;
    ScanDir = ScanDir { 1 } ( 3 ) ;
    DirCheck = strfind ( ScanDir , 'Down' ) ;

    % Rotate all data according to whether up or down scan.
    if isempty ( DirCheck { 1 } ) == 1 % i.e. if it is an up scan.

        for j = 1 : size ( DataArray , 2 )
            % Make all data 2D.
            DataArray { 1 , j , i } = reshape ( DataArray { 1 , j , i } , SamplesPerLine , NumLines ) ;
            % Over-rotate up scans by 270 degrees and mirror them so they start at the top.
            DataArray { 1 , j , i } = flip ( rot90 ( DataArray { 1 , j , i } , 3 ) , 2 );

        end

    else % i.e if it is a down scan.

        for j = 1 : size ( DataArray , 2 )
            % Make all data 2D.
            DataArray { 1 , j , i } = reshape ( DataArray { 1 , j , i } , SamplesPerLine , NumLines ) ;
            % Rotate down scans by 90 degrees.
            DataArray { 1 , j , i } = rot90 ( DataArray { 1 , j , i } ) ;

        end

    end

    %% Section to combine trace/retrace data.

    % Interpolate/concatenate trace and retrace current and voltage into a single
    % column, sorting by sets of length SamplesPerLine, and generate time
    % values. Assumes first current column is trace data and second is
    % retrace, because Nanoscpe appears to label all exported data as
    % trace, so this cannot be determined from the data headers.
    
    k = 0 ; % Index for user to check correct import direction.
    TraceFirst = 1 ; % 0 if retrace data is first. Trace first as default.

    while k < 1

        for j = 1 : NumLines

            clear CurrentLine
            clear VoltageLine

            % Remove starting artefact. Unsure of cause, could be an export
            % issue.
            if size ( iV , 2 ) == 2 && j == 1 && mode ( DataArray { 1 , iV ( 2 ) , i } ( 1 , : ) ) < 0.0019

                DataArray { 1 , iC ( 2 ) , i } ( 1 , : ) = 0 ;
                DataArray { 1 , iV ( 2 ) , i } ( 1 , : ) = 0 ;

            end

            % Combine trace and retrace data lines. Retrace data is generally
            % first, although sometimes trace is first. Not sure why yet.

            if TraceFirst == 0

                CurrentLine ( : , 2 ) = horzcat ( flip (...
                    DataArray { 1 , iC ( 2 ) , i } ( j , : ) , 2 ) ,...
                    DataArray { 1 , iC ( 1 ) , i } ( j , : ) )' ;

                elseif TraceFirst == 1

                CurrentLine ( : , 2 ) = horzcat (...
                    DataArray { 1 , iC ( 1 ) , i } ( j , : ) ,...
                    flip ( DataArray { 1 , iC ( 2 ) , i } ( j , : ) , 2 ) )' ;

            end

            if size ( iV , 2 ) == 2 % Check whether trace and retrace voltages are captured.

                if TraceFirst == 0

                VoltageLine ( : , 1 ) = abs ( horzcat ( flip (...
                    DataArray { 1 , iV ( 2 ) , i } ( j , : ) , 2 ) ,...
                    DataArray { 1 , iV ( 1 ) , i } ( j , : ) ) )' ;

                elseif TraceFirst == 1

                VoltageLine ( : , 1 ) = abs ( horzcat (...
                    DataArray { 1 , iV ( 1 ) , i } ( j , : ) ,...
                    flip ( DataArray { 1 , iV ( 2 ) , i } ( j , : ) , 2 ) ) )' ;

                end

            else

                % Just duplicate voltage data if only a single direction was
                % captured.

                VoltageLine ( : , 1 ) = abs ( horzcat (...
                    DataArray { 1 , iV , i } ( j , : ) ,...
                    ( DataArray { 1 , iV , i } ( j , : ) .* 0 ) +...
                    DataArray { 1 , iV , i } ( j , end ) ) )' ;

            end

            % Generate sampling times and time axis.
            StartTime = CurrentAndTime ( end , 1 ) ;
            FinishTime = StartTime + ( 2 * SamplesPerLine * SampleTime ( i , 1 ) ) - SampleTime ;
            CurrentLine ( : , 1 ) = ( StartTime : SampleTime ( i , 1 ) : FinishTime )' ;
            % Include voltage recorded in file.
            VoltageLine ( : , 2 ) = ApplV ;

            % Concatenate new values to previous.
            CurrentAndTime = vertcat ( CurrentAndTime , CurrentLine ) ;
            VoltageAndTime = vertcat ( VoltageAndTime , VoltageLine ) ;
            CurrentPerFile { i } = vertcat ( CurrentPerFile { i } , CurrentLine ) ;
            VoltagePerFile { i } = vertcat ( VoltagePerFile { i } , VoltageLine ) ;

            if j == 1

                SampleTime ( i , 2 ) = CurrentLine ( 1 , 1 ) ; % Start time of file.

            elseif j == NumLines

                SampleTime ( i , 3 ) = CurrentLine ( end , 1 ) ; % End time of file.

            end

        end
        
        figure ;
        plot ( VoltageAndTime ( : , 1 ) ) ; % Quickly look at data to check.
        set ( gcf , 'Color' , 'w' , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;
        DirAnswer = questdlg ( 'Is data import direction correct?' , 'Check import direction.' ) ;
        close ;
        switch DirAnswer
            
            case 'Yes'
                
                k = 1 ; % End loop and proceed.
                
            case 'No'
                
                % Clear previous data and try again with opposite import direction.
                CurrentAndTime = zeros ( 1 , 2 ) ;
                VoltageAndTime = zeros ( 1 , 2 ) ;
                TraceFirst = abs ( TraceFirst - 1 ) ;
                
            case 'Cancel'
                
                return ;
                
        end
                            
    end

    %% Remove starting values.

    CurrentPerFile { i } ( 1 , : ) = [] ;
    VoltagePerFile { i } ( 1 , : ) = [] ;

    PrevTime = FinishTime ;

end

%% Remove starting values.
CurrentAndTime ( 1 , : ) = [] ;
VoltageAndTime ( 1 , : ) = [] ;

%% Produce smoothed current and voltage data.
FilterWidth = 256 ;
CurrentAndTime ( : , 3 ) = smooth ( CurrentAndTime ( : , 2 ) , FilterWidth , 'sgolay' , 1 ) ;
VoltageAndTime ( : , 3 ) = smooth ( VoltageAndTime ( : , 1 ) , FilterWidth , 'sgolay' , 1 ) ;

%% Additional metrics; ressitance and conductance.

% ResistanceAndTime = VoltageAndTime ;
% ResistanceAndTime ( : , 1 ) = ResistanceAndTime ( : , 1 ) ./ ( CurrentAndTime ( : , 2 ) .* 1E-9 ) ;
% ResistanceAndTime ( : , 2 ) = ResistanceAndTime ( : , 2 ) ./ ( CurrentAndTime ( : , 2 ) .* 1E-9 ) ;
% ConductanceAndTime = ResistanceAndTime ;
% G0 = 2 * ( 1.60217662E-19 ^ 2 ) / 6.62607004E-34 ;
% ConductanceAndTime ( : , 1 ) = ( 1 ./ ( ConductanceAndTime ( : , 1 ) ) ) ./ G0 ;
% ConductanceAndTime ( : , 2 ) = ( 1 ./ ( ConductanceAndTime ( : , 2 ) ) ) ./ G0 ;
% SmoothResistance = smooth ( ResistanceAndTime ( : , 2 ) , 32 , 'sgolay' , 1 ) ;

%% Determine whether measurement is constant current or constant voltage.

FeedbackCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'TUNAList' ) ) ) ) ) ;
FeedbackType = textscan ( CalArray { FeedbackCell ( 1 ) } , '%s' ) ;
FeedbackType = FeedbackType { 1 } { 4 } ;

SampleThickness = 1.1 ; % Thickness in 10s of nm, to convert readily to MV/cm.

if  ~contains ( FeedbackType , 'Dis' ) == 1 % For constant current mode, i.e. TUNAList is Enabled.

    % Constant current sensitivity.
    ConstantSensCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'Sens. tunaSens' ) ) ) ) ) ;
    ConstantSens = textscan ( CalArray { ConstantSensCell ( 1 ) } , '%s' ) ;
    ConstantSens = str2double ( ConstantSens { 1 } { 4 } ) ;
    % Current setpoint.
    SetPointCell = find ( not ( cellfun ( 'isempty' , ( strfind ( CalArray ,...
        'TunaFeedbackSetpoint' ) ) ) ) ) ;
    SetPoint = textscan ( CalArray { SetPointCell ( 1 ) } , '%s' ) ;
    SetPoint = sprintf ( '%g' , ConstantSens * str2double ( SetPoint { 1 } { 7 } ) ) ;
    DataTitle = strcat ( 'Constant' , { ' ' } , SetPoint , { ' ' } , 'nA stress' ) ;
    CharArray ( 1 ) = abs ( str2double ( SetPoint ) ) ;

    figure ;

    yyaxis left
    semilogy ( CurrentAndTime ( : , 1 ) , CurrentAndTime ( : , 3 ) , '-' ) ;
    hold on
    ylabel ( 'Current/nA' ) ;
    xlabel ( 'Time/s' ) ;

    yyaxis right
    plot ( CurrentAndTime ( : , 1 ) , VoltageAndTime ( : , 1 ) ./ SampleThickness , 'LineWidth' , 2 ) ;
    % Voltage converted to field.
    hold on
    ylabel ( 'Electric Field/MVcm^{-1}' ) ;
    set ( gca , 'FontSize' , 24 ) ;
    set ( gcf , 'Color' , 'w' , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

    title ( DataTitle ) ;

    DataOut = horzcat ( CurrentAndTime , VoltageAndTime ) ;
    save ( strcat ( DataTitle { 1 } , '.txt' ) , 'DataOut' , '-ascii' ) ;
    saveas ( gcf , strcat ( DataTitle { 1 } , '.tif' ) ) ;

else % For constant voltage mode, i.e. TUNAList is Disabled.

    %% Section to find and note onset of TDDB if present. This will be done
    % automatically, though it may just give unusable information for non-TDDB
    % measurements.

    %% Indices for voltage on and breakdown times to use.
    tV = 2 ;
    tI = 1 ; % Setting to 3 seems to avoid artefacts and catch some part of the first BD peak.
    OnsetCurrent = 1 ; % Define onset current in nA.

    DataTitle = strcat ( num2str ( round ( max ( VoltageAndTime ( : , 1 ) ) , 2 ) ) , ' V applied' ) ;
    disp ( strcat ( 'At a constant' , { ' ' } , DataTitle ) ) ;
    CharArray ( 1 ) = round ( max ( VoltageAndTime ( : , 1 ) ) , 2 ) ;

    % Determine current and time values relative to when voltage swithced on.
    % Display the time at which the current threshold is exceeded.

    if size ( iV , 2 ) == 1
        % Best pracrise to use current spike artifact if only one voltage
        % direction is recorded, because the actual onset of voltage might
        % be lost in the un-recorded direction.
        IStart = findchangepts ( VoltageAndTime ( : , 1 ) , 'MaxNumChanges' , 20 ) ;
        VStart = findchangepts ( CurrentAndTime ( IStart ( tV )- 1024 :...
            IStart ( tV ) + 1024 , 2 ) , 'MaxNumChanges' , 20 ) + IStart ( tV ) - 1024 ;

    else
        % Use change in voltage if both directions are recorded.
        VStart = findchangepts ( VoltageAndTime ( : , 1 ) , 'MaxNumChanges' , 20 ) ;

    end

    VIndicator = [ CurrentAndTime( VStart( tV ) , 1 ) , min( VoltageAndTime ( : , 1 ) ) ;...
        CurrentAndTime( VStart ( tV ) , 1 ) , max( VoltageAndTime ( : , 1 ) ) ] ;

    TDDBI = find ( CurrentAndTime ( VStart ( tV ) + 1 : end , 2 ) >= OnsetCurrent ) ;

    %% Section for plotting V or I vs time and saving plots.

    figure ;
    subplot ( 2 , 1 , 1 ) ;

    yyaxis left
    semilogy ( CurrentAndTime ( : , 1 ) , CurrentAndTime ( : , 2 ) , '-' ) ;
    hold on
    semilogy ( CurrentAndTime ( : , 1 ) , CurrentAndTime ( : , 3 ) , 'y-' ) ;
    ylabel ( 'Current/nA' ) ;
    xlabel ( 'Time/s' ) ;

    yyaxis right
    % Voltage converted to field.
    plot ( CurrentAndTime ( : , 1 ) , VoltageAndTime ( : , 1 ) ./ SampleThickness , 'LineWidth' , 2 ) ;
    hold on
    plot ( VIndicator ( : , 1 ) , VIndicator ( : , 2 ) , 'b--' ) ;
    ylabel ( 'Electric Field/MVcm^{-1}' ) ;

    if isempty ( TDDBI ) == 0

        IIndicator = [ CurrentAndTime( VStart ( tV ) + TDDBI( tI ) , 1 ) ,...
            min( CurrentAndTime ( : , 2 ) ) + 0.01 ; CurrentAndTime( VStart ( tV )...
            + TDDBI ( tI ) , 1 ) , max( CurrentAndTime ( : , 2 ) ) ] ;

        yyaxis left
        semilogy ( IIndicator ( : , 1 ) , IIndicator ( : , 2 ) , 'r--' ) ;

        OnsetTime = CurrentAndTime ( VStart ( tV ) + TDDBI ( tI ) , 1 ) -...
        CurrentAndTime ( VStart ( tV ) , 1 ) ;
        disp ( strcat ( 'Time to breakdown:' , { ' ' } ,...
        num2str ( OnsetTime ) , { ' ' } , 's' ) ) ;
        title ( strcat ( DataTitle , { ', ' } , 'Time to breakdown:' , { ' ' } ,...
        num2str ( OnsetTime ) , { ' ' } , 's' ) ) ;

    else

        disp ( 'No breakdown.' ) ;
        DataTitle = strcat ( DataTitle , { ', no breakdown' } ) ;
        DataTitle = DataTitle { 1 } ;
        title ( DataTitle ) ;

        TDDBI = size ( CurrentAndTime , 1 ) - VStart ( tV ) - 2048 ; % Arbitrary value so that SILC may subsequently still be fitted.

    end

    set ( gca , 'FontSize' , 24 ) ;
    set ( gcf , 'Color' , 'w' , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

    DataOut = horzcat ( CurrentAndTime , VoltageAndTime ) ;
    save ( strcat ( DataTitle , '.txt' ) , 'DataOut' , '-ascii' ) ;

end

%% Section for peak analysis of signal or part(s) of signal.

figure ;
subplot ( 3 , 1 , 1 ) ;
yyaxis left
semilogy ( CurrentAndTime ( : , 1 ) , CurrentAndTime ( : , 2 ) ) ;
ylabel ( 'Current/nA' ) ;
xlabel ( 'Time/s' ) ;
yyaxis right
plot ( CurrentAndTime ( : , 1 ) , VoltageAndTime ( : , 1 ) , 'LineWidth' , 2 ) ;
ylabel ( 'Voltage/V' ) ;
set ( gca , 'FontSize' , 24 ) ;
set ( gcf , 'Color' , 'w' , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

BoundDefault = { '0' , '0' , '0' , '0.2' , '0' , '0.1' } ;
a = 0 ; % Peak analysis index.
while a < 1 % Choose appropriate crop/sampling region and perform peak analysis.

    % Get user input of time to sample between.
    BoundPrompt = { 'Enter lower bound time:' , 'Enter upper bound time:' ,...
            'Enter voltage peak height threshold (in V above local background):' ,...
            'Enter voltage peak prominence threshold (in V):' ,...
            'Enter current peak height threshold (in nA above local background):' ,...
            'Enter current peak prominence threshold (in nA):' } ;
    BoundTitle = 'Sampling region input.' ;
    BoundAnswer = inputdlg ( BoundPrompt , BoundTitle , 1 , BoundDefault ) ;
    LowerTime = str2double ( BoundAnswer { 1 } ) ;
    UpperTime = str2double ( BoundAnswer { 2 } ) ;
    VPeakThresh = str2double ( BoundAnswer { 3 } ) ;
    VPromThresh = str2double ( BoundAnswer { 4 } ) ;
    CPeakThresh = str2double ( BoundAnswer { 5 } ) ;
    CPromThresh = str2double ( BoundAnswer { 6 } ) ;

    [ LMin , LowerBound ] = min ( abs ( CurrentAndTime ( : , 1 ) - LowerTime ) ) ;
    [ UMin , UpperBound ] = min ( abs ( CurrentAndTime ( : , 1 ) - UpperTime ) ) ;

    BoundDefault = BoundAnswer' ;

    % Get part of signal for analysis and plot.
    clear AnalysisSection ;
    AnalysisSection ( : , 1 ) = CurrentAndTime ( LowerBound : UpperBound , 1 ) ;
    AnalysisSection ( : , 2 ) = CurrentAndTime ( LowerBound : UpperBound , 2 ) ;
    AnalysisSection ( : , 3 ) = VoltageAndTime ( LowerBound : UpperBound , 1 ) ;
    subplot ( 3 , 1 , 1 ) ;
    yyaxis left
    plot ( AnalysisSection ( : ,1 ) , AnalysisSection ( : , 2 ) , 'LineWidth' , 2 ) ;
    ylabel ( 'Current/nA' ) ;
    xlabel ( 'Time/s' ) ;
    yyaxis right
    plot ( AnalysisSection ( : , 1 ) , AnalysisSection ( : , 3 ) );
    ylabel ( 'Voltage/V' ) ;
    set ( gca , 'FontSize' , 24 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    title ( DataTitle ) ;

    % Get and check frequency information.
    [ FMin , FLower ] = min ( abs ( SampleTime ( : , 3 ) - LowerTime ) ) ;
    [ FMax , FUpper ] = min ( abs ( SampleTime ( : , 2 ) - UpperTime ) ) ;

    FreqLower = SampleTime ( FLower , 1 ) ;
    FreqUpper = SampleTime ( FUpper , 1 ) ;

    % To notify user if the sample rate varies through the measurement.
    if FreqLower ~= FreqUpper

        disp 'Sample rates are inconsistent across analysis region.' ;

    end

    subplot ( 3 , 1 , 2 ) ;
    findpeaks ( AnalysisSection ( : , 3 ) , ( 1 / FreqLower ) , ...
        'MinPeakProminence' , VPromThresh , 'Threshold' , VPeakThresh ) ;
    ylabel ( 'Voltage/V' ) ;
    xlabel ( 'Time/s' ) ;
    set ( gca , 'FontSize' , 24 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    title ( 'Voltage peaks' ) ;

    subplot ( 3 , 1 , 3 ) ;
    findpeaks ( AnalysisSection ( : , 2 ) , ( 1 / FreqLower ) , ...
        'MinPeakProminence' , CPromThresh , 'Threshold' , CPeakThresh ) ;
    ylabel ( 'Current/nA' ) ;
    xlabel ( 'Time/s' ) ;
    set ( gca , 'FontSize' , 24 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    title ( 'Current peaks' ) ;

    hold off

    % Let user choose whether to continue working on FFTs.
    PeakChoice = questdlg ( 'Continue sampling or complete process:' ,...
        'Sampling options:' , 'Repeat' , 'Stop session' , 'Complete' ,...
        'Repeat' ) ;

    switch PeakChoice

        case 'Repeat'

        case 'Stop session'

            h = 1 ;
            a = 1 ;
            OutputData = { ' ' } ; % Output something to avoid error.
            CharArray = zeros ( 1 , 4 ) ;

        case 'Complete'

            a = 1 ;

            saveas ( gcf , strcat ( 'PeaksOut' , '.tif' ) ) ;

            % Once analysis section and desired promience threshold chosen,
            % generate a range of thresholds about the chosen value.
            VIterationRange = 0 : VPromThresh / 4 : 5 * VPromThresh ;
            OutputData = zeros ( numel ( VIterationRange ) , 9 ) ;
            SetThreshold = VPromThresh ; % Retain set prominence threshold.

            for i = 1 : numel ( VIterationRange )

                % To iterate through a range of prominence thresholds.
                VPromThresh = VIterationRange ( i ) ;

                [ VoltagePeaks , VoltageLocs , VoltageWidths , VoltageProms ] = ...
                    findpeaks ( AnalysisSection ( : , 3 ) , ( 1 / FreqLower ) , ...
                    'MinPeakProminence' , VPromThresh , 'Threshold' , VPeakThresh ,...
                    'WidthReference' , 'HalfProm' ) ;

                %% Get histogram data for spread/frequency of peaks/spikes.

                % Spike separation.
                VPeakSeps = diff ( VoltageLocs ) ;
                if numel ( VPeakSeps ) < 2
                    
                    OutputData ( i , 1 ) = 0 ;
                    OutputData ( i , 2 ) = 0 ;
                    OutputData ( i , 3 ) = 0 ;
                
                else
                    
                    SepDist = fitdist ( VPeakSeps , 'Gamma' ) ;
                    VPeakSpread = pdf ( SepDist , VPeakSeps ) ;
                    OutputData ( i , 1 ) =  max ( VPeakSpread ) ;
                    OutputData ( i , 2 ) = SepDist.a ;
                    OutputData ( i , 3 ) = SepDist.b ;
                
                end
                
                % Spike widths.
                if numel ( VoltageWidths ) < 2

                    OutputData ( i , 4 ) = 0 ;
                    OutputData ( i , 5 ) = 0 ;
                    OutputData ( i , 6 ) = 0 ;
                
                else
                    
                    WidthDist =  fitdist ( VoltageWidths , 'Gamma' )  ;
                    VWidthSpread = pdf ( WidthDist , VoltageWidths ) ;
                    OutputData ( i , 4 ) = max ( VWidthSpread ) ;
                    OutputData ( i , 5 ) = WidthDist.a ;
                    OutputData ( i , 6 ) = WidthDist.b ;

                end
                
                % Spike prominences.
                if numel ( VoltageProms ) < 2
                    
                    OutputData ( i , 7 ) = 0 ;
                    OutputData ( i , 8 ) = 0 ;
                    OutputData ( i , 9 ) = 0 ;
                    
                else
                    
                    PromsDist =  fitdist ( VoltageProms , 'Gamma' )  ;
                    VPromsSpread = pdf ( PromsDist , VoltageProms ) ;
                    OutputData ( i , 7 ) = max ( VPromsSpread ) ;
                    OutputData ( i , 8 ) = PromsDist.a ;
                    OutputData ( i , 9 ) = PromsDist.b ;

                end
                
                % Spike heights.
                if numel ( VoltagePeaks ) < 2
                    
                    OutputData ( i , 10 ) = 0 ;
                    OutputData ( i , 11 ) = 0 ;
                    OutputData ( i , 12 ) = 0 ;
                    
                else
                    
                    HeightsDist =  fitdist ( VoltagePeaks , 'Gamma' )  ;
                    VHeightsSpread = pdf ( HeightsDist , VoltagePeaks ) ;
                    OutputData ( i , 10 ) = max ( VHeightsSpread ) ;
                    OutputData ( i , 11 ) = HeightsDist.a ;
                    OutputData ( i , 12 ) = HeightsDist.b ;

                end
                
            end

    end

end

%% Section for FFT analysis of signal or part(s) of signal.

% Index whether or not to perform FFT analysis.
% if  strfind ( FeedbackType , 'Dis' )  > 0 % For constant voltage mode, i.e. TUNAList is Enabled.
%     
%     f = 1 ;
% 
% else
%     
%     f = 0 ;
% 
%     figure ;
%     subplot ( 2 , 1 , 1 ) ;
%     yyaxis left
%     semilogy ( CurrentAndTime ( : , 1 ) , CurrentAndTime ( : , 2 ) ) ;
%     ylabel ( 'Current/nA' ) ;
%     xlabel ( 'Time/s' ) ;
%     yyaxis right
%     plot ( CurrentAndTime ( : , 1 ) , VoltageAndTime ( : , 1 ) , 'LineWidth' , 2 ) ;
%     ylabel ( 'Voltage/V' ) ;
%     set ( gca , 'FontSize' , 24 ) ;
%     set ( gcf , 'Color' , 'w' , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;
%     
% end
% 
% BoundDefault = { '0' , '0' , '5' , '50' } ;
% 
% while f < 1 % Choose appropriate crop/sampling region and perform FFT.
%     
%     % Get user input of time to sample between.
%     BoundPrompt = { 'Enter lower bound time:' , 'Enter upper bound time:' ,...
%             'Enter sampling window 1 in seconds:' ,...
%             'Enter sampling window 2 in seconds:' } ;
%     BoundTitle = 'Sampling region input.' ;
%     BoundAnswer = inputdlg ( BoundPrompt , BoundTitle , 1 , BoundDefault ) ;
%     LowerTime = str2double ( BoundAnswer { 1 } ) ;
%     UpperTime = str2double ( BoundAnswer { 2 } ) ;
%     FFTWindow1 = str2double ( BoundAnswer { 3 } ) ;
%     FFTWindow2 = str2double ( BoundAnswer { 4 } ) ;
%     
%     [ LMin , LowerBound ] = min ( abs ( CurrentAndTime ( : , 1 ) - LowerTime ) ) ;
%     [ UMin , UpperBound ] = min ( abs ( CurrentAndTime ( : , 1 ) - UpperTime ) ) ;
%     
% 	BoundDefault = BoundAnswer' ;
% 
%     % Ensure the sampling window has an even sample length.
%     if rem ( UpperBound - LowerBound , 2 ) ~= 0
%         
%         UpperBound = UpperBound + 1 ;
%         
%     end
%     
%     clear AnalysisSection ;
%     AnalysisSection ( : , 1 ) = CurrentAndTime ( LowerBound : UpperBound , 1 ) ;
%     AnalysisSection ( : , 2 ) = VoltageAndTime ( LowerBound : UpperBound , 1 ) ;
%     subplot ( 2 , 1 , 1 ) ;
%     yyaxis left
%     plot ( CurrentAndTime ( LowerBound : UpperBound , 1 ) ,...
%         CurrentAndTime ( LowerBound : UpperBound , 2 ) , 'LineWidth' , 2 ) ;
%     ylabel ( 'Current/nA' ) ;
%     xlabel ( 'Time/s' ) ;
%     yyaxis right
%     plot ( AnalysisSection ( : , 1 ) , AnalysisSection ( : , 2 ) );
%     ylabel ( 'Voltage/V' ) ;
%     set ( gca , 'FontSize' , 24 ) ;
%     set ( gcf , 'Color' , 'w' ) ;
%     
%     [ FMin , FLower ] = min ( abs ( SampleTime ( : , 3 ) - LowerTime ) ) ;
%     [ FMax , FUpper ] = min ( abs ( SampleTime ( : , 2 ) - UpperTime ) ) ;
%     
%     FreqLower = SampleTime ( FLower , 1 ) ;
%     FreqUpper = SampleTime ( FUpper , 1 ) ;
%     
%     % To notify user if the sample rate varies through the measurement.
%     if FreqLower ~= FreqUpper
%         
%         disp 'Sample rates are inconsistent across analysis region.' ;
%         
%     end
%     
%     % Check the signal has an even length and append the final value if
%     % not.
%     FullSignalLength = size ( AnalysisSection , 1 ) ;
%     EvenCheck = rem  ( FullSignalLength , 2 ) ;
%     if EvenCheck ~= 0
%         
%         AnalysisSection = vertcat ( AnalysisSection ,...
%             AnalysisSection ( end , : ) ) ;
%         FullSignalLength = size ( AnalysisSection , 1 ) ;
%         
%     end
%     
%     % Perform full FFT and get into frequency domain.
%     FullFFT = fft ( AnalysisSection ( : , 2 ) ) ; % FFT of current in analysis section.
%     % Transform parameters.
%     FullFrequency = ( 1 ./ FreqLower ) * ( 0 : ( FullSignalLength / 2 ) ) / FullSignalLength ;
%     FullSpec2 = abs ( FullFFT / FullSignalLength ) ; % Two-sided spectrum.
%     FullSpec1 = FullSpec2 ( 1 : ( FullSignalLength / 2 ) + 1 ) ; % Single-sided spectrum.
%     FullSpec1 ( 2 : end - 1 ) = 2 .* FullSpec1 ( 2 : end - 1 ) ;
%             
%     % Determine parameters to chop up signal into smaller chunks.
%     % Number of samples in each sampling window.
%     WindowSamples1 = round ( FFTWindow1 / FreqLower ) ;
%     WindowSamples2 = round ( FFTWindow2 / FreqLower ) ;
%     
%     % Ensure window sampling lengths are even.
%     EvenCheck = rem  ( WindowSamples1 , 2 ) ;
%     if EvenCheck ~= 0
%         
%         WindowSamples1 = WindowSamples1 + 1 ;
%         
%     end
%     
%     EvenCheck = rem  ( WindowSamples2 , 2 ) ;
% 
%     if EvenCheck ~= 0
%         
%         WindowSamples2 = WindowSamples2 + 1 ;
%         
%     end
%     
%     % Output vectors.
%     WindowOutput1 = zeros ( ( WindowSamples1 / 2 ) + 1 , 1 ) ;
%     WindowOutput2 = zeros ( ( WindowSamples2 / 2 ) + 1 , 1 ) ;
%     % Transform parameters.
%     Frequency1 = ( 1 ./ FreqLower ) * ( 0 : ( WindowSamples1 / 2 ) ) / WindowSamples1 ;
%     Frequency2 = ( 1 ./ FreqLower ) * ( 0 : ( WindowSamples2 / 2 ) ) / WindowSamples2 ;
%     
%     c = 0 ; % Chopping index.
%     c1 = 0 ; % Window 1 counting index.
%     c2 = 0 ; % Window 2 counting index.
%     w1 = 1 ; % Window 1 index. Windex?
%     w2 = 1 ; % Window 2 windex. Windex!
%     p1 = 1 ; % Count number of window 1s.
%     p2 = 1 ; % Count number of window 2s.
%     
%     while c < 1
% 
%         % Check that next iteration won't tery to sample beyond the length
%         % of the data for first sampling window.
%         if w1 + WindowSamples1 -1 <= size ( AnalysisSection , 1 ) && c1 < 1
%             
%             % Capture sampling window.
%             AnalysisWindow1 = AnalysisSection ( w1 : w1 + WindowSamples1 - 1 , : ) ;
%             % Perform  FFT on sampling windows and get into frequency domain.
%             FFTSection1 = fft ( AnalysisWindow1 ( : , 2 ) ) ; % FFT of current in analysis section.
%             Window1Spec2 = abs ( FFTSection1 / WindowSamples1 ) ; % Two-sided spectrum.
%             Window1Spec1 = Window1Spec2 ( 1 : ( WindowSamples1 / 2 ) + 1 ) ; % Single-sided spectrum.
%             Window1Spec1 ( 2 : end - 1 ) = 2 .* Window1Spec1 ( 2 : end - 1 ) ;
% 
%             % Append data to output vector.
%             WindowOutput1 = WindowOutput1 + Window1Spec1 ;
%             w1 = w1 + WindowSamples1 ;
%             p1 = p1 + 1 ;
% 
%         else
%             
%             c1 = 1 ; % Sampling window 1 has completed running through the data.
%         
%         end
%         
%         % Perform for second sampling window.
%         if w2 + WindowSamples2 - 1 <= size ( AnalysisSection , 1 ) && c2 < 1
% 
%             AnalysisWindow2 = AnalysisSection ( w2 : w2 + WindowSamples2 - 1 , : ) ;
%             FFTSection2 = fft ( AnalysisWindow2 ( : , 2 ) ) ;
%             Window2Spec2 = abs ( FFTSection2 / WindowSamples2 ) ;
%             Window2Spec1 = Window2Spec2 ( 1 : ( WindowSamples2 / 2 ) + 1 ) ;
%             Window2Spec1 ( 2 : end - 1 ) = 2 .* Window2Spec1 ( 2 : end - 1 ) ;
% 
%             WindowOutput2 = WindowOutput2 + Window2Spec1 ;
%             w2 = w2 + WindowSamples2 ;
%             p2 = p2 + 1 ;
%         
%         else
%             
%              c2 = 1 ;
%         
%         end
%         
%         % End if both windows have traversed the full data.
%         if c1 == 1 && c2 == 1
% 
%             c = 1 ;
% 
%         end
% 
%     end    
% 
%     % Scale data according to number of sampling iterations.
%     WindowOutput1 = WindowOutput1 ./ p1 ;
%     WindowOutput2 = WindowOutput2 ./ p2 ;
%     
%     subplot ( 2 , 1 , 2 ) ;
%     semilogx ( FullFrequency , FullSpec1 ) ;
%     ylabel ( 'Amplitude/a.u.' ) ;
%     xlabel ( 'Frequency/Hz' ) ;
%     hold on
%     semilogx ( Frequency1 , WindowOutput1 ) ;
%     semilogx ( Frequency2 , WindowOutput2 ) ;
%     legend ( 'Full data' , char ( strcat ( BoundAnswer { 3 } , { 's window' }  ) ) ,...
%         char ( strcat ( BoundAnswer { 4 } , { 's window' } ) ) ) ;
%     set ( gca , 'FontSize' , 24 ) ;
%     set ( gcf , 'Color' , 'w' ) ;
%     
%     title ( DataTitle ) ;
%     hold off
%     
%     % Let user choose whether to continue working on FFTs.
%     FFTChoice = questdlg ( 'Continue sampling or complete process:' ,...
%         'Sampling options:' , 'Repeat' , 'Repeat with hold' , 'Complete' ,...
%         'Repeat' ) ;
%     
%     switch FFTChoice
%         
%         case 'Repeat'
%             
%         case 'Repeat with hold'
%             
%             h = 1 ;
%             hold on
%             
%         case 'Complete'
%             
%             f = 1 ;
%             
%             saveas ( gcf , strcat ( 'FFTOut' , '.tif' ) ) ;
%
%     end
%     
% end

% Record characteristics of measurement.
CharArray ( 2 ) = max ( CurrentAndTime ( : , 2 ) ) ;
CharArray ( 3 ) = max ( VoltageAndTime ( : , 1 ) ) ;
CharArray ( 4 ) = mean ( VoltageAndTime ( LowerBound : UpperBound , 1 ) ) ;

end

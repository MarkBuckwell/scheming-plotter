function [ NormData ] = XPSNormaliser ( FullDataArray )

%% Normalises ARXPS depth profile data to a distint photoemission peak, asuming
% that more than one peak has been sampled in a measurement.

NameArray = FullDataArray ( 1 , : ) ;

for i = 1 : length ( NameArray )
   
    NameArray{ i } = NameArray { i } ( 1 ) ;
    
end

% Determines which photoemission should be normalised to.
NormIndex = listdlg ( 'PromptString' , 'Choose data to normalise to:' , 'ListString' , NameArray ) ;

N = numel ( OutputHeaders ) ;
NFC = length ( DataArray ) ;
NDepth = NFC / NAngles ; % Number of depths sampled.
BackSubArray = repmat ( { '' } , 1 , N - 2 , NFC ) ;
% Array for background-subtracted data, with 1 column removed as this is
% the background subtracted. Energy scale is maintained in order to produce
% correctly-scaled integral areas. Configuration might require adjustment
% for different output files. Assumes first column is binding energy data.

Background = find ( ~cellfun ( @isempty , ( strfind ( OutputHeaders , 'Back' ) ) ) ) ;
% Finds column index of background data to subtract. Errors will occur if
% two headers contain the string 'Back', so this might need modify,
% depending on the input data.

OutputHeaders ( Background ) = [] ;
OutputHeaders ( 1 ) = [] ;
NOut = numel ( OutputHeaders ) ;

for i = 1 : NFC
    
    l = 2 ; % Index for subtraction.
            % Assumes second column is first intensity data for background
            % subtraction.
    
    for k = 1 : N - 1

        if k == Background - 1
            
            % Do nothing in the case that the indexed column is the
            % backtround data.
            
        else
    
        BackSubArray{ 1 , l - 1 , i } = DataArray { 1 , k + 1 , i } - DataArray { 1 , Background , i } ;
        l = l + 1 ; % Only incrememnt subtraction index when indexed column in not background data.
        
        end
        
    end
    
end

BackSubMatrix = zeros ( 1 , NOut , NFC ) ; % Matrix for integrated areas.

for i = 1 : NFC
   
    for k = 1 : N - 2 
        
        BackSubMatrix ( 1 , k , i ) = trapz ( flipud ( DataArray { 1 , 1 , i } ) , flipud ( BackSubArray { 1 , k , i } ) ) ;
        % Data flipped, assuming that binding energy scale is reversed
        % (otherwise integrals come out negative).
        
    end
    
end

%% Processes data to account for configuration of sampled angles and depths.

ARDataArray = repmat ( { '' } , NAngles , NDepth ) ;
PlotData = zeros ( NAngles , NDepth , NOut ) ;
DataTransfer = zeros ( NAngles , NDepth ) ;

for i = 1 : NOut
   
    for k = 1 : NFC
    
        ARDataArray{ k } = BackSubMatrix ( 1 , : , k ) ;
        DataTransfer ( k )  = ARDataArray { k } ( i ) ;
    end
    
    PlotData ( : , : , i ) = DataTransfer ;
    
end

end

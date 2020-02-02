clear
clc

%% Imports and processes XPS data columns in order to extract area ratios etc.
% Used to plot trends in peak components as a function of depth and
% emission angle. Can be used on ASCII data exported from CasaXPS.

addpath ( cd ) ;
DataSetup = inputdlg ( { 'Number of datasets (distinct regions) to import:' ,...
    'Number of angles sampled:' , 'Smallest angle:' , 'Angle increment:' } ,...
    'Set up data import.' , [ 1 35 ] , { '3' , '12' , '22.5 ' , '5' } ) ;

FullDataArray = repmat ( { '' } , 2 , str2double ( DataSetup { 1 } ) ) ;

NameArray = repmat ( { '' } , 1 , str2double ( DataSetup { 1 } ) ) ;

Angles = str2double ( DataSetup { 3 } ) : str2double ( DataSetup { 4 } )...
    : str2double ( DataSetup { 3 } ) + ( str2double ( DataSetup { 4 } )...
    * ( str2double ( DataSetup { 2 } ) - 1 ) ) ; % Photoemission angles sampled.

for i = 1 : str2double ( DataSetup { 1 } )
    
    % Import data from files corresponding to a distinct photoemission peak (i.e. a region spectrum).
    [ DataArray , OutputHeaders , DataName ] = XPSFileSelector ( ) ;

    % Gets absolute peak areas by subtracting background from data and fitting envelope.
    [ PlotData , OutputHeaders ] = XPSSubtractor ( DataArray , OutputHeaders , str2double ( DataSetup { 2 } ) ) ;

    % Puts data and headers into full array.
    OutputHeaders{ 1 } = DataName ;
    NameArray{ i } = DataName ;
    FullDataArray{ 1 , i } = OutputHeaders ;
    FullDataArray{ 2 , i } = PlotData ;
    
end

if str2double ( DataSetup { 1 } ) == 1

    Depths = 0 : 1 : ( size ( PlotData , 2 ) ) - 1 ;
    figure ( 1 ) ;
    hold on
    PlotColumns = round ( length ( OutputHeaders ) / 2 ) ;
    PlotRows = round ( length ( OutputHeaders ) / PlotColumns ) ;

    for j = 1 : numel ( OutputHeaders )

        subplot ( PlotRows , PlotColumns , j ) ;
        imagesc ( Angles , fliplr ( Depths ) , rot90 ( PlotData ( : , : , j ) ) ) ;
        title ( OutputHeaders { j } ) ;
        ylabel ( 'Etch time/s' ) ;
        xlabel ( 'Emission angle/degrees' ) ;
        colorbar ;

    end
    
else
    
    [ NormData ] = XPSNormaliser ( FullDataArray ) ;
    
end

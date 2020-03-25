clear
clc

%% Imports constant voltage or current data, running through processing,
% to ouptut data of current and voltage in time, as well as peak detection
% parameters (distibution of separations, prominences, widths, heights).
% Input data should be 2D CAFM data (Bruker/Veeco) wherein the scan size is
% 0, producing instead data in time, rather than space. Script will
% continue looping until the user cancels importing an additional file.

% DataTitles = [ { 'Setpoint/nA or V' } , { 'Peak current/nA' } , { 'Peak voltage/V' } ,...
%     { 'Mean settled voltage/V' } ] ;
% DataTitles = [ { 'Modal spike separation/s' } , { 'Gamma shape' } , { 'Gamms scale' } ,...
%     { 'Modal spike width/s' } , { 'Gamma shape' } , { 'Gamma scale' } ,...
%     { 'Modal spike prominence/V' } , { 'Gamma shape' } , { 'Gamma scale' } ,...
%     { 'Modal spike height/V' } , { 'Gamma shape' } , { 'Gamma scale' } ] ;
Characteristics = zeros ( 1 , 4 ) ;
i = 0 ; % To continue adding files/data.
j = 0 ;
NFC = 1 ; % File index.
while i < 1

    [ OutputData , CharArray , SetThreshold , VIterationRange , h ] =...
        CAFMTimeVsCurrentFunction ( ) ; % Get data.
    % Append measurement data.
    FullDataMetrics { NFC } = OutputData ;
    Characteristics = vertcat ( Characteristics , CharArray ) ;
    NFC = NFC + 1 ; % Increment file index.
    
    if h == 1 % End session.
        i = 1 ;
    end
    close all

end

NFC = NFC - 1 ; % Correct number of files.
Characteristics ( 1 , : ) = [ ] ; % Delete first, unneeded row.

%% Sort data into ascending order by setpoint.
[ SortedChars , CharArrange ] = sortrows ( Characteristics ) ;
% Condutance values in nS.
SortedChars ( : , 5 ) = SortedChars ( : , 1 ) ./ SortedChars ( : , 4 ) ;
% Array for sorted data.
SortedData = zeros ( 26 , 17 , NFC ) ;
for i = 1 : NFC
    
    SortedData ( : , : , CharArrange ( i ) ) = FullDataMetrics { i } ;

end

%% Plot basic voltage and current characteristics.
set ( 0 , 'DefaultFigureWindowStyle' , 'Docked' ) ;

figure ;
loglog ( SortedChars ( : , 1 ) , SortedChars ( : , 2 ) , 'o' ) ;
xlabel ( 'Current bias/nA' ) ;
ylabel ( 'Peak current/nA' ) ;
set ( gca , 'FontSize' , 14 ) ;
set ( gcf , 'Color' , 'w' ) ;

figure ;
semilogx ( SortedChars ( : , 1 ) , SortedChars ( : , 3 ) , 'o' ) ;
xlabel ( 'Current bias/nA' ) ;
ylabel ( 'Peak voltage/V' ) ;
set ( gca , 'FontSize' , 14 ) ;
set ( gcf , 'Color' , 'w' ) ;

figure ;
semilogx ( SortedChars ( : , 1 ) , SortedChars ( : , 4 ) , 'o' ) ;
xlabel ( 'Current bias/nA' ) ;
ylabel ( 'Mean settled voltage/V' ) ;
set ( gca , 'FontSize' , 14 ) ;
set ( gcf , 'Color' , 'w' ) ;

figure ;
loglog ( SortedChars ( : , 1 ) , SortedChars ( : , 5 ) , 'o' )
xlabel ( 'Current bias/nA' ) ;
ylabel ( 'Conductance/nS' ) ;
set ( gca , 'FontSize' , 14 ) ;
set ( gcf , 'Color' , 'w' ) ;

%% Plot distribution metrics.
% Get row corresponding to set promience threshold.

figure ;
for i = 1 : numel ( VIterationRange )
    
    loglog ( SortedChars ( : , 1 ) ,...
         ( reshape ( SortedData ( i , 2 , : ) , NFC , 1 ) ) ./...
         reshape ( SortedData ( i , 3 , : ) , NFC , 1 ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal time between spikes/s' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    
end

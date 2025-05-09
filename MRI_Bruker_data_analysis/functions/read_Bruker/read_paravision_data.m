function scan = read_paravision_data(path_2dseq)
% based on the read_paravision_data_PV6 function by Tim Klasen, 2013
% tiny modifications by Samo Lasic, 2022

%%  DATEI-AUSWAHL
%
%   Wird kein Funktionsargument geliefert, so öffnet sich ein Fenster, um
%   eine Datei auszuwählen
%
%   Es werden drei Pfade als String angelegt:   - Hauptpfad der Messung
%                                               - Nebenpfad (Ordner der
%                                                 Reko)
%                                               - Ort der 2dseq-Datei
%
%   Wird ein ungültiger Pfad übergeben, wird die Funktion mit einer Fehler-
%   meldung abgebrochen
%
%%

% file separator
fs = filesep;


if nargin == 0
    [filename, scan.info.paths.directory_sub] = uigetfile(...
        {'*.*', 'File (*.*)';},...
        '2dSeq-File auswählen!');

    if strcmp(filename,'2dseq') == 1 %|| strcmp(filename,'d3proc') == 1
        scan.info.paths.directory_2dseq = [scan.info.paths.directory_sub filename];
        bs = strfind(scan.info.paths.directory_sub,fs);
        disp((length(bs)));
        scan.info.paths.directory_main = scan.info.paths.directory_sub; %(1:bs(length(bs)+ 0));

    else
        error('Input was no 2dseq-file! Please enter a valid path!')
    end
else
    bs = strfind(path_2dseq,fs);
    if strcmp(path_2dseq(bs(end)+1:end),'2dseq') == 1 %|| strcmp(path_2dseq(bs(end)+1:end),'d3proc') == 1
        scan.info.paths.directory_sub = path_2dseq(1:end-5);
        scan.info.paths.directory_2dseq = path_2dseq;
        bs = strfind(scan.info.paths.directory_sub,fs);
        scan.info.paths.directory_main = scan.info.paths.directory_sub(1:bs(length(bs)-2));
    else
        error('Input was no 2dseq-file! Please enter a valid path!')
    end
end

%%  LISTE DER EINZULESENDEN DATEIEN
%
%   in "files-main" werden die Dateien angegeben, die in dem Haupt-
%   verzeichnis stehen (scan.info.paths.directory_main)
%
%   in "files_sub" werden die Dateien angegeben, die in dem Reko-
%   verzeichnis stehen (scan.info.paths.directory_sub)
%
files_main = {...
    'acqp';...
    'method';...
    };

files_sub = {...
    'd3proc';...
    'id';...
    'procs';...
    'reco';...
    'visu_pars';...
    };

files = [files_main; files_sub];

%%  EINLESEN DER BRUKER-DATEIEN
%
%   Pro Durchlauf der äußeren For-Schleife wird genau eine Datei
%   ausgelesen, sofern diese existiert

for index = 1:length(files)

    %   Datei aus dem Hauptpfad wird geöffnet und gelesen (sofern diese
    %   existiert)
    if index <= length(files_main)
        scan.info.paths.(files{index}) = [scan.info.paths.directory_main files{index}];
        scan.info.fid.(files{index}) = fopen(scan.info.paths.(files{index}),'r');

        if (scan.info.fid.(files{index}) == -1)
            warning('MATLAB:FileOpen', [files{index} '-file could not be open or does not exist!'])
        else
            temp = fread(scan.info.fid.(files{index}));
            scan.info.strings.(files{index}) = char(temp');
            fclose(scan.info.fid.(files{index}));
            act_row = scan.info.strings.(files{index});
        end

        %   Datei aus dem Reko-Pfad wird geöffnet und gelesen (sofern diese
        %   existiert)
    else
        scan.info.paths.(files{index}) = [scan.info.paths.directory_sub files{index}];
        scan.info.fid.(files{index}) = fopen(scan.info.paths.(files{index}),'r');

        if (scan.info.fid.(files{index}) == -1)
            warning('MATLAB:FileOpen', [files{index} '-file could not be open or does not exist!'])
        else
            temp = fread(scan.info.fid.(files{index}));
            scan.info.strings.(files{index}) = char(temp');
            fclose(scan.info.fid.(files{index}));
            act_row = scan.info.strings.(files{index});
        end

    end


    if (scan.info.fid.(files{index}) ~= -1)
        % RR erzeugt eine Matrix mit den Position der Zeichen "##" UND "##$"
        % RRD erzeugt eine Matrix mit den Postionen des Zeichen "##$"
        RR = strfind(scan.info.strings.(files{index}),'##');
        RRD = strfind(scan.info.strings.(files{index}),'##$');

        RRDi = 1;
        RRi = 1;

        line_break = sprintf('\n');

        while ~isempty(RRD)



            if length(RRD) == 1

                scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) = char(act_row(RRD(1):RR(2)-1));
                scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                    strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), line_break, '');
                scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                    strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), '##$', '');

                DD = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'$$');

                if ~isempty(DD)
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) = ...
                        scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(1:DD-1);
                end

                % entfernt Klammerausdrücke nach dem "="-zeichen
                left_parenthesis = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'(');
                right_parenthesis = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),')');
                length_row = length(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]));

                if ~isempty(left_parenthesis) && length_row ~= right_parenthesis(1)
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                        [scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(1:left_parenthesis(1)-1)...
                        scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(right_parenthesis(1)+1:end)];
                end

                % entfernt spitze Klammern
                scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                    strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), '<', '');
                scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                    strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), '>', '');

                % sucht das "="-Zeichen und unterteilt den String in Parameter und den Wert
                equal_sign = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'=');
                parameter = scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(1:equal_sign);
                parameter_value = scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(equal_sign+1:end);

                if ~isempty(str2num(scan.info.header.(files{index}).(['RRD' num2str(RRDi)])...
                        (strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'=')+1:end)))
                    eval(['scan.(files{index}).' parameter '[' parameter_value '];'])
                else
                    %                         eval(['scan.(files{index}).' parameter '''' parameter_value ''';'])
                end




                break

            else

                %% Die nächste Zeile fängt mit ""##" oder "##$" an
                %  Schleife formatiert die Zeile, so dass am Ende ein Wert
                %  ausgegeben wird

                % "##$"-Zeichen wurde gefunden
                if RR(1) == RRD(1)

                    % erstellt eine einzelne Zeile und entfernt Zeilenumbrüche/##$-Zeichen
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) = char(act_row(RRD(1):RRD(2)-1));
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                        strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), line_break, '');
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                        strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), '##$', '');

                    % entfernt die Zeichen ab $$
                    DD = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'$$');
                    if ~isempty(DD)
                        scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) = ...
                            scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(1:DD-1);
                    end

                    % entfernt Klammerausdrücke nach dem "="-zeichen
                    left_parenthesis = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'(');
                    right_parenthesis = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),')');
                    length_row = length(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]));

                    if ~isempty(left_parenthesis) && length_row ~= right_parenthesis(1)
                        scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                            [scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(1:left_parenthesis(1)-1)...
                            scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(right_parenthesis(1)+1:end)];
                    end

                    % entfernt spitze Klammern
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                        strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), '<', '');
                    scan.info.header.(files{index}).(['RRD' num2str(RRDi)]) =...
                        strrep(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]), '>', '');

                    % sucht das "="-Zeichen und unterteilt den String in Parameter und den Wert
                    equal_sign = strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'=');
                    parameter = scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(1:equal_sign);
                    parameter_value = scan.info.header.(files{index}).(['RRD' num2str(RRDi)])(equal_sign+1:end);

                    if ~isempty(str2num(scan.info.header.(files{index}).(['RRD' num2str(RRDi)])...
                            (strfind(scan.info.header.(files{index}).(['RRD' num2str(RRDi)]),'=')+1:end)))
                        eval(['scan.(files{index}).' parameter '[' parameter_value '];'])
                    else
                        eval(['scan.(files{index}).' parameter '''' parameter_value ''';'])
                    end


                    RRDi = RRDi + 1;
                    act_row = act_row(RRD(2):end);


                    % "##"-Zeichen wurde gefunden
                else
                    RRi = RRi + 1;
                    act_row = act_row(RR(2):end);
                end

            end

            RR = strfind( act_row,'##');
            RRD = strfind( act_row,'##$');



        end
    end

end

%% Einlesen des Bilddatensatzes
if scan.info.fid.d3proc == -1
    fid=fopen(scan.info.paths.directory_2dseq,'r','l');
    nSlices = scan.visu_pars.VisuCoreFrameCount;
    %    nSlices = scan.visu_pars.VisuCoreSlicePacksSlices;

    % a fix for DWS data
    if isfield(scan.visu_pars,'VisuAcqPhaseEncSteps')

        % a temporary fix for scanner-applied drift correction
        if prod(scan.visu_pars.VisuAcqSize == scan.reco.RECO_ft_size) == 0
            phases = scan.reco.RECO_ft_size(2);
            reads = scan.reco.RECO_ft_size(1);
        else
            phases = scan.visu_pars.VisuAcqPhaseEncSteps;
            reads  = scan.visu_pars.VisuAcqSize(1);
        end

    else
        phases = 1;
        reads  = scan.visu_pars.VisuAcqSize(1);
    end



    scan.dataset = zeros([reads , phases, nSlices]);
    size(scan.dataset);
    if strcmp(scan.visu_pars.VisuCoreWordType, '_16BIT_SGN_INT')
        encoding = 'int16';
    elseif strcmp(scan.visu_pars.VisuCoreWordType, 'ip_short')
        encoding = 'short';
    elseif strcmp(scan.visu_pars.VisuCoreWordType, 'ip_int')
        encoding = 'int';
    end
    for slice = 1:nSlices
        for phase = 1:phases
            a = fread(fid, reads, encoding);
            scan.dataset(1:reads, phase, slice) = a;
        end
    end
    fclose(fid);
end

scan.dataset = permute(scan.dataset,[3 2 1]);
scan.datasetOrder = 'Slice x Phase x Read';


end
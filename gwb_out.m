function [C] = gwb_out(data,datatype,param)
%C = GWB_OUT(DATA,DATATYPE,PARAM) parses text output files produced by
%Geochemist's Work Bench based on user input
%
%   DATA = Data from GWB output file stored as cell array with 1 line per
%          row
%
%   DATATYPE = section of file where output paramater of interest should be
%              located, entered as a text string. Only one datatype allowed
%              for each function call.
%
%   PARAM = Exact name of GWB output paramater to look for stored as a
%           string in a cell array. Search string must match format in
%           output file (e.g., for sulfate, must use the form SO4-- NOT
%           SO42-), but is not case sensitive. Multiple PARAM values
%           allowed for each datatype, separated by a comma.
%   
%   TOT option = One additional feature allowed for aqueous speices is the
%                calculation of the total amount of a species; this
%                option is invoked by adding 'TOT' in front of the species
%                name. Note this method will not work for independently
%                totalling different oxidation states of the same species
%   --> Need to add units option so user can specify different options for
%       units for a given parameter
%   --> Just handling this by adding new cases as needed

%% FILE MANAGEMENT BLOCK
%--> Now handled in executive function to avoid redundant operations
while 0
    %Use only for manual file management
    %Read in entire output file
    filename = ['React_output_kinetic_knallgas_2nd_order_9.30e-01.txt'];
    fid = fopen(filename);
    output = textscan(fid,'%s','delimiter','\n','whitespace','');
    fclose(fid);
    %Extract contents of cell to form suitable for GETLINES
    data = output{1};
end

%% HEADERS BLOCK

%If a data is empty, end the script
if isempty(data)
    error('No data. Likely cause is output file not found.')
end

%Create empty array for function return if error encountered
C = {};

%Get line number vectors for section headers of interest. 
%
%NOTES
% 
% -Search string is not case sensitive but must otherwise match string in
%  gwb file.
% -Use line numbers returned here to operate on individual chunks of file.
% -Create format specifiers based on datatype. 
% -Format identifiers must allow for parameter name/value pairs and should
%  be constructed so only one output goes into each output cell

switch lower(datatype)
    case 'system'
        top = getlines(data,'Step #');
        %NOTE: nernst block may or may not be present.
        %      Several possible workarounds:
        
        %1. Sort ascending on line numbers such that moving through chunks
        %   in the order they are printed out. The line number of the chunk
        %   printed first is always in the first column, and this is what
        %   is returned for each row when only one indexed is used as in
        %   the for-loop over the 'top' vector below.
        b1 = getlines(data,'Nernst redox');
        b2 = getlines(data,'Reactants');
        bot = sort([b1 b2],2); %works with one or more non-empty vectors
        
        %2. Universal fix for all output sections would be to add print
        %   commands in GWB template so all sections are always present
        %   (might not work because some sections such as Nernst and Redox
        %   don't seem to be recognized GWB print categories):
        %       print species long
        %       print surfaces long
        %       print saturations long
        %       print gases long
        %       print basis long
        %       print orig_basis long
        %       print elements long
        %       print reactions long
        
        %3. Search for the last entry in the system block (Water type, then
        %   add one to its line number. This method untested; not sure if
        %   it's reproducible for all output files
        %bot = getlines(data,'Water type') + 1;
        
        for i = 1:length(param)
            %Choose format specifier based on paramater
            %--> Make param indexing uniform, thouth it doesn't seem to
            %    matter whether it's written with {} or ()
            if strcmpi(param{i},'Step #')
                %Get Step # label/value
                fs{i} = '%s %*s %f %*s %*s %*f';
            elseif strcmpi(param{i},'Time')
                %Get Time label/value (in seconds)
                fs{i} = '%s %*s %f %*s %*s %*s';
                %NOTE: Could also get days, but would need to use this
                % format specifier:
                %      fs = '%*s %*s %*f %*s %f %*s';
                % and add an exception to the Extract-and-store
                % section with these options for the textscan statement:
                %      'Delimiter',{' ','(',')'},'MultipleDelimsAsOne',1
            elseif strcmpi(param{i},'Temperature')
                %Get temperature label/value
                fs{i} = '%s %*s %f %*s %*s %*s %*f %*s';
            elseif strcmpi(param{i},'Pressure')
                %Get pressure label/value
                fs{i} = '%*s %*s %*f %*s %s %*s %f %*s';
            elseif strcmpi(param(i),'pH')
                %Get pH label/value pair
                fs{i} = '%s %*s %f %*s %*s %*s %*f';
            elseif strcmp(param(i),'Solvent mass')
                %Get solvent mass label/value pair
                fs{i} = '%s %*s %*s %f %*s';
            elseif strcmp(param(i),'Solution mass')
                %Get solvent mass label/value pair
                fs{i} = '%s %*s %*s %f %*s';
            elseif strcmp(param(i),'Solution density')
                %Get solution density label/value pair
                fs{i} = '%*s %s %*s %f %*s';
            elseif strcmp(param(i),'Fluid volume')
                %Get solution volume label/value pair
                fs{i} = '%*s %s %*s %f %*s';
            else
                fprintf(1,'Need new format specifier for %s\n',param{i});
                return
            end
        end
    case 'nernst'
        top = getlines(data,'Nernst redox');
        bot = getlines(data,'Reactants');
        fprintf(1,'Need new format specifier for %s\n',datatype);
    case 'reactant'
        top = getlines(data,'Reactants');
        %NOTE: kinetic block may or may not be present
        %NOTE: mineral section header changes depending on whether minerals
        %      are isolated from system
        b1 = getlines(data,'Kinetic');
        b2 = getlines(data,'Minerals in system');
        b3 = getlines(data,'Minerals isolated');
        %See 'system' case for explanation of sorting approach
        bot = sort([b1 b2 b3],2); %works with one or more non-empty vectors
        fs{1} = '%s %*f %*f %f %*f';
    case 'minerals'
        %NOTE: mineral section header changes depending on whether minerals
        %      are isolated from system
        t1 = getlines(data,'Minerals in system');
        t2 = getlines(data,'Minerals isolated');
        top = sort([t1;t2]);   %works only if all but one vectors are empty
        bot = getlines(data,'Aqueous species');
        %Get mineral name and moles isolated from system
        fs{1} = '%s %f %*f %*f %*f';
        %Get mineral name and volume of mineral isolated from system
        %fs{1} = '%s %*f %*f %*f %f';
    case 'aqueous'
        %NOTE: using 'aqueous' instead of 'activity' for legacy code
        %--> add a case to get activity coefficients. Can't just extract
        %    this value at the same time as concentration b/c program is
        %    designed to return only a single number per chemical species
        %    and modifying it could wreak havoc in unforeseen ways.
        top = getlines(data,'Aqueous species');
        bot = getlines(data,'Mineral saturation states');
        %Get species name and log activity
        fs{1} = '%s %*f %*f %*f %f';
    case 'molality'
        %Same section of output as 'aqueous', but gets molality instead of
        %log activity
        top = getlines(data,'Aqueous species');
        bot = getlines(data,'Mineral saturation states');
        %Get species name and molality
        fs{1} = '%s %f %*f %*f %*f';
    case 'saturation'
        top = getlines(data,'Mineral saturation states');
        bot = getlines(data,'Gases');
        fs{1} =  '%s %f %s %f';
    case 'gases'
        top = getlines(data,'Gases');
        bot = getlines(data,'Original basis');
        fprintf(1,'Need new format specifier for %s\n',datatype);
    case 'current'
        top = getlines(data,'Basis components');
        %NOTE: If present, current basis may be followed by original basis
        %      or elemental composition depending on print options set in
        %      reaction script
        b1 = getlines(data,'Original basis');
        b2 = getlines(data,'Elemental composition');
        %See 'system' case for explanation of sorting approach
        bot = sort([b1 b2],2); %works with one or more non-empty vectors
        %Get species name and total moles in fluid
        %NOTE: this should be the same as total moles if precipitation is
        %      turned off 
        fs{1} = '%s %*f %f %*f';
    case 'original'
        top = getlines(data,'Original basis');
        bot = getlines(data,'Elemental composition');
        %Get species name and total moles in fluid
        %NOTE: this should be the same as total moles if precipitation is
        %      turned off 
        fs{1} = '%s %*f %f %*f';
    case 'element'
        top = getlines(data,'Elemental composition');
        bot = getlines(data,'Step #');
        %Get rid of first entry (top of file) & add entry for end of file
        bot(1) = [];
        bot(length(bot)+1) = length(data);
        fprintf(1,'Need new format specifier for %s\n',datatype);
    otherwise
        error(['Unrecognized data type: ''%s''. Verify that datatype '...
               'is correct and that requested and bracketing data '...
               'chunks are set to print to output file'],datatype);
end

if ~exist('fs','var')
    error('Format Specifier Missing')
end

%--> ERROR CHECKING FOR DIFFERENT NUMBER OF SECTION HEADERS
%--> Currently handling this on case by case basis using SORT

%% OUTPUT PARSING

%Loop through all entries in PARAM vector to get line #'s for each
for i=1:length(param)
    %Set some initial switches
    tot_conc = false;
    
    %Store data label
    C{1,i} = param{i};
    
    %Check to see if the total amount of an aqueous species is requested
    if strcmp(datatype,'aqueous') && length(param{i}) > 3
        tot_test = param{i}(1:3);
        if strcmpi('TOT',tot_test)
            %Set switch to true; strip out the 'tot' designator
            %--> Need to account for weak acids such as H2S and H2CO3 which
            %    will have different names when dissociated
            %--> Consider just using Basis components section
            %    which store the total for any basis species!!!
            tot_conc = true;
            param{i} = param{i}(4:length(param{i}));
        end
    end     
    
    %Loop through all data blocks in the file
    for j=1:length(top)
        %Set format identifier
        if strcmp(datatype,'system')
            %Each system parameter requires separate format identifier
            fslocal = fs{i};
        else
            fslocal = fs{1};
        end
        
        %Set start/end indices for data chunk
        if strcmp(datatype,'system')
            d1 = top(j);  
        else
            d1 = top(j)+1;  %Did this b/c mg/kg was getting caught
                            %in search for Mg
        end
        d2 = bot(j);
        chunk = data(d1:d2);
        
        %If dealing with saturation indices, scrub the lines of any text
        %information about saturation state
        if strcmp(datatype,'saturation')
            chunk = regexprep(chunk,'s/sat| sat ','');
        end
        
        %Get index within chunk for line containing parameter
        di = getlines(chunk,param{i});
        
        %If more than one match found and not looking for total
        %concentration, find exact match
        if length(di)>1 && ~tot_conc
            for k=1:length(di)
                lab_val = textscan(chunk{di(k)},fslocal);
                if strcmp(datatype,'saturation')
                    %Have to check both name columns for exact match
                    if strcmp(lab_val{1},param{i}) ||...
                       strcmp(lab_val{3},param{i})
                        di = di(k);
                        break
                    end                    
                else
                    if strcmp(lab_val{1},param{i})
                        di = di(k);
                        break
                    end
                end
            end
        end
        
        %Extract and store quantity of interest
        if isempty(di) || (length(di)>1 && ~tot_conc)
            %If matching algorithm did not produce expected result
            %Uncomment warning for debug purposes
            %warning(sprintf(['Did not find exact match for '...
            %    'string ''%s'' in Chunk #: %d'],param{i},j))
            quantity(j,1) = NaN;
        elseif length(di)>1 && tot_conc
            %For getting total concentrations
            quantity(j,1) = 0;
            %Change format specifier to molalities, which can be added
            %VERIFIED: get the same number whether adding molalities or
            %          dividing original basis moles by solvent mass,
            %          though second option is probably easier
            fslocal = '%s %f %*f %*f %*f';
            for k=1:length(di)
                lab_val = textscan(chunk{di(k)},fslocal);
                quantity(j,1) = quantity(j,1) + lab_val{2};
            end
            %***Adding molalities, which is ok, but need to MULTIPLY by
            %total mass of solvent to get total moles of basis species***
            %TOT option would require calling gwb_out recursively to get
            %mass of solution for current reaction step
            %--> Not possible with current program. Instead, get conc and
            %    mass vectors separately, then multiply the two together.
            %--> OR divide vector of original basis moles by total mass of
            %    solution to get molality total basis species at each step.
            %    Should also work for each species involved in redox pair
            %    when that redox pair is coupled.
        elseif strcmp(datatype,'saturation')
            %For mineral saturation index
            lab_val = textscan(chunk{di},fslocal);
            %Remove empty cells, accounts for the fact that there are two
            %minerals per line and the desired mineral is not always in the
            %same spot.
            lab_val = lab_val(~cellfun('isempty',lab_val));
            lvi = find(cellfun(@(x) strcmpi(x,param{i}), lab_val));
            strMatch = lab_val{lvi}{1};
            quantity(j,1) = lab_val{lvi+1};
        else
            %For everything else
            lab_val = textscan(chunk{di},fslocal);
            strMatch = lab_val{1}{1};
            quantity(j,1) = lab_val{2};
        end
        
        %Do final check to make sure the right parameter was gotten. Not
        %applicable for 'system' datatypes and tot_conc option because
        %these are allowed to have inexact matches with the search string
        if length(di)==1 && ~tot_conc && ~strcmpi(datatype,'system') && ...
                ~strcmpi(param{i},strMatch)
            %Uncomment warning for debug purposes
            %warning(sprintf(['Found string ''%s'' does not match search '...
            %   'string ''%s'' in Chunk #: %d'],strMatch,param{i},j))
            quantity(j,1) = NaN;
            clear strMatch  %Just to be safe, may not be necessary
        end
        
        %--> ERROR CHECKING FOR:
        %   - DATATYPE/PARAMETER MISMATCH
        %   - BAD FORMAT SPECIFIER
        
    end
    %Store output
    C{2,i} = quantity;
    %--> Consider returning data as Table:
    %
    %First need to convert parameter names (which must match strings in
    %GWB output files) to suitable variable names
    % colNames = matlab.lang.makeValidName(C(1,:));
    %
    %Create Table
    % T = table(C{2,:},'VariableNames',colNames);
    %This accessed using dot notation and meaningful variable names, e.g.:
    %T.SolutionMass for Solution Mass
    %
    %-->Downside to this is in order to loop through a table with an
    %   unknown number of columns, would need to use indices into
    %   VariableNames(): T.(T.Properties.VariableNames{1})
    %   This is no more clear than the cell array approach in use now.
    %   Also, other codes use cell functions to deal with output.
end

%End of function GWB_OUT
end
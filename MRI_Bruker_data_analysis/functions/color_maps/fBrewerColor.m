function [col N] = fBrewerColor(num,flip)
    % num is sequential number 
    % set optional flip='flip' argument to reverse the colors
    
    % Brewer Color Maps
    % selected for color blind, printing, photocopy, LCD friendly
    % http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
    %
    % Color Blind Friendly: This icon indicates that a given color scheme will not confuse people
    % with red-green color blindness. Red-green color blindness affects approximately 8 percent
    % of men and 0.4 percent of women, although its severity varies and so some schemes will have a "?"
    % indicating it may be a problem for some, but not all folks with color vision impairment.
    %
    % Color Printing Friendly: Suitable for desktop color printing (based on our color printers -
    % your mileage will vary and you should check on your own printer). We checked Matchprint proofs
    % for all of these schemes and they are all printed in a 2003 journal article and the book
    %     Designing Better Maps by Cynthia Brewer. CMYK specs are as close to press-ready as is reasonable.
    %
    % Photocopy Friendly: This indicates that a given color scheme will withstand black and white
    % photocopying. Diverging schemes can not be photocopied successfully.
    % Differences in lightness should be preserved with sequential schemes.
    %
    % Laptop (LCD) Friendly: This icon indicates that a given color scheme is suitable for viewing
    % on a laptop LCD display. Small, portable LCD monitors tend to wash-out colors which results in
    % noticeable differences from computer-to-computer.
    
    
    %super safe
    BC{1} = [254,230,206;253,174,107;230,85,13];
    BC{end+1} = [222,235,247; 158,202,225;49,130,189];
    BC{end+1} = [240,240,240; 189,189,189;99,99,99];
    BC{end+1} = [241,163,64; 247,247,247; 153,142,195];
    
    BC{end+1} = [254,240,217;253,204,138;252,141,89;215,48,31];
   

    
    % safe for color blind, but not for black and white printing
    BC{end+1} = [254,237,222; 253,190,133; 253,141,60; 217,71,1];
    BC{end+1} = [239,243,255; 189,215,231; 107,174,214; 33,113,181];
    BC{end+1} = [247,247,247; 204,204,204; 150,150,150; 82,82,82];
    
    BC{end+1} = [254,237,222; 253,190,133; 253,141,60;230,85,13; 166,54,3];
    BC{end+1} = [239,243,255; 189,215,231; 107,174,214; 49,130,189; 8,81,156];
    BC{end+1} = [247,247,247; 204,204,204; 150,150,150; 99,99,99; 37,37,37];
    
    BC{end+1} = [254,237,222; 253,208,162; 253,174,107; 253,141,60; 230,85,13; 166,54,3];
    BC{end+1} = [239,243,255;198,219,239;158,202,225;107,174,214;49,130,189;8,81,156];
    BC{end+1} = [247,247,247;217,217,217;189,189,189;150,150,150;99,99,99;37,37,37];
    
    BC{end+1} = [254,237,222;253,208,162;253,174,107;253,141,60;241,105,19;217,72,1;140,45,4];
    BC{end+1} = [239,243,255;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,69,148];
    BC{end+1} = [247,247,247;217,217,217;189,189,189;150,150,150;115,115,115;82,82,82;37,37,37];
    
    BC{end+1} = [255,245,235;254,230,206;253,208,162;253,174,107;253,141,60;241,105,19;217,72,1;140,45,4];
    BC{end+1} = [247,251,255;222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,69,148];
    BC{end+1} = [255,255,255;240,240,240;217,217,217;189,189,189;150,150,150;115,115,115;82,82,82;37,37,37];
    
    BC{end+1} = [255,245,235;254,230,206;253,208,162;253,174,107;253,141,60;241,105,19;217,72,1;166,54,3;127,39,4];
    BC{end+1} = [247,251,255;222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107];
    BC{end+1} = [255,255,255;240,240,240;217,217,217;189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0];
    
    col = BC{num}/255;
    
    if nargout > 1
        
        N = size(col,1);
    end
    
    if (exist('flip', 'var'))
        if strcmp(flip,'flip') 
        col = flipud(col);
        end
        
    end
    
end










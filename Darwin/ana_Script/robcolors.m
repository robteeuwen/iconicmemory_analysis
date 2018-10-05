function M = robcolors(varargin)

% this function breeds fancy colors; 
% 
% input(colormap name, number of colors)
% - moonrisekingdom (default) (3 colors default)  

if isempty(varargin)
    m = 'moonrisekingdom';
else
    if strcmp(varargin{1},'moonrisekingdom')
        m = 'moonrisekingdom';
    end
end

colnum = 0; 

if length(varargin) > 1
    colnum = varargin{2}; 
end

if strcmp(m,'moonrisekingdom')
    T = [133, 175,   163
         225,   176, 36 
         142,   25,   30
         133, 175,   163]./255; 
    Ts = [0 64 191 255];
    
    if ~colnum
        colnum = 3;
    end
    
    if colnum > size(T,1)
        M = interp1(Ts/255,T,linspace(0,1,colnum));
    else
        M = T(1:colnum,:); 
    end
end



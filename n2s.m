function [str, ixe] = n2s(str, ixe, val, prec, sfx, right2left)
% Description
%
% codegen-able num2str. accepts scalars only!
%
% str[1*N char]                 - char array buffer
% ixe[1*1 int32]                -  index of last used character (0 for initial)
% val[1*1 double]               -  value to be printed
% prec[1*1 double/single]       -  decimal precision (positive for fractions, 0 for integers, negative for trailing zeros up to 10^0 order)
% sfx[1*N char]                 - string to suffix the written data
%
% modifies str buffer and returns index of last active character in buffer
%
%
% TODO: maybe use this for integers only, and the rest write in scientific
% (X.YeZ). This will remove the need of varying precision (extreme low/exteme high values will be written/read just fine)

%%
coder.varsize('str', [1, 1e4], [0 1]);

% single value to buffer string *uses precision, and adds space in end
if nargin<6
    right2left = false;
end
if coder.target('matlab') && isinteger(val)
    val = double(val);
end
sg          = sign(val);
val         = val*sg;
if val>1e16||val<1e-16
    % cpp protection against exploding values (=gets corrupted on high values)
    val = 0;
end
len         = int32(0);
lenmax      = length(str);
% length of buffer for single (scalar) value (since dynamic range of double is ~16 orders magnitude, there should be no need to display values above 17 characters.)
% **THIS function dhould not be used for values above(or below) 1e16/1e-16 (for that scientific, and not simple decimal, notation should be used)
str1        = char(zeros(1, 21)); % set max #characters to 21 in any case
isvb40      = false;%prec>=0; % value before prec=0 reached (decimal .)
isfirst     = true;
% protection against stray values from epsilon (max possible precision for given value)
prec        = double(prec);
prec        = min(prec, -1-floor(log10(eps(val)))); 
if prec<0
    % should be separate from while val>0 || prec>=0;end
    % added as protection from eps(val) interferance with LSB values in extreme dynamic ranges, = 
    % while>end loop writes wrong(random) numbers
    for i=1:-prec
        len             = len + int32(1);
        str1(len)       = 48;
    end
%     warning('MATLAB:CustomNum2Str:ExtremeDynamicRange', 'requesting n2s with value above possible dynamic range of double. Replaced out-of-range values with trailing zeros')
end
val             = round(val/10^-prec);
while val>0 || prec>=0
    
    rm              = mod(val, 10); % extract lsb value
    isvb40          = isvb40 || rm>0; % flag for non-zero values (ignores trailing zeros in decimal fractions)
    if prec==0 && isvb40 && ~ isfirst
        len             = len + int32(1);
        str1(len)       = '.';
        isfirst         = false;
    end
    if isvb40||prec<=0 % write @ prec<=0 OR prec>0 after a non-zero LSB
        len             = len + int32(1);
        str1(len)       = 48 + uint8(rm);
        isfirst         = false;
    end
    prec            = prec-1;
    val             = (val - rm)/10;
end
if sg<0
    len             = len + int32(1);
    str1(len)       = '-';
end
% if nargin>4
%     L = length(sfx);
%     for i=1:L
%         str(ixe+i)      	= sfx(i);
%     end
%     ixe            	= ixe+L;
% end


rem2write = int32(0); % remaining characters till buffer end
if right2left
    rem2write = ixe;
else
    rem2write = lenmax - ixe;
end
if rem2write<len
    % write error into string buffer
    str(ixe) = 'E';
    return
else
    if right2left
        rem2write = ixe;
        str(ixe+1:-1:ixe-len+2) = str1(1:len);
    else
        rem2write = lenmax - ixe;
        str(ixe+[1:len])= str1(len:-1:1);
    end
end

ixe            	= ixe+len;
if ~right2left && nargin>4 %suffix will work only on left-to-right writing
    L = length(sfx);
    for i=1:L
        str(ixe+i)      	= sfx(i);
    end
    ixe            	= ixe+L;
end

end


%% test validity

%{
N =  1000;


n = randi(14, N, 1);
v = rand(N, 1);

val = v.*10.^n;

str = char(zeros(1, 100e3));
ixe = 0;
for i=1:N
    [str, ixe]=n2s_exp(str, ixe, val(i), 16);
end

str2 = char(zeros(1, 100e3));
ixe2 = 0;
for i=1:N
    [str2, ixe2]=n2s(str2, ixe2, val(i), 16);
end


str3 = sprintf('%1.16f ', val);

% scientific vs decimal precision
(str2num(str)-str2num(str2))./val';

% n2s vs sprintf [both decimal precision]
(str2num(str3)-str2num(str2))./val';

%}
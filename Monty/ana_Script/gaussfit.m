function [y,params] = gaussfit(x,r,Starting,hardbase)

%fminsearch approach
%Fits a basic GAussian

%x is input value (time, space etc)
%r = response;

%Starting Parameters are [mean,std,Gain and baseline]
%No further normalisation is required.

%If hardbase is inputted then no baseline shift is allowed and the baseline
%is set to hardbase

options=optimset('Display','off');
 
if nargin<4
    %Baseline fitted
   
    params=fminsearch(@gaussfitinline,Starting,options,x,r);
    %Remake best fitting Guass
    mx = params(1);
    sx = params(2);
    G = params(3);
    b = params(4);
    
    y = G.*normpdf(x,mx,sx)+b;
else
    %Basleine restricted
    params=fminsearch(@gaussfitinlinenobase,Starting(1:3),options,x,r,hardbase);
    %Remake best fitting Guass
    mx = params(1);
    sx = params(2);
    G = params(3);
    
    y = G.*normpdf(x,mx,sx)+hardbase;
end


return

function sse = gaussfitinline(params,x,r)

mx = params(1);
sx = params(2);
G = params(3);
b = params(4);

m = G.*normpdf(x,mx,sx)+b;

sse = sum(sum((m-r).^2));

return

function sse = gaussfitinlinenobase(params,x,r,hardbase)

mx = params(1);
sx = params(2);
G = params(3);

m = G.*normpdf(x,mx,sx)+hardbase;

sse = sum(sum((m-r).^2));

return
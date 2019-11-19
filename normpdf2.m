function y=normpdf2(x1,x2,mu,cm)
% Bivariate normal PDF
%
% J. Ditterich, 10/01
%
% y = normpdf2 (x1,x2,mu,cm)
%
% y is a matrix (of the same size as X1 and X2) containing the values of the PDF.
%
% x1 defines the locations on the first axis. It can be a matrix of arbitrary shape.
% x2 defines the locations on the second axis. It must have the same size as X1.
% mu defines the location of the maximum. It must be a 2-element vector.
% cm is the covariance matrix.

% History:
% released on 10/4/01 as part of Toolbox V 0.6 Beta
% additional checks added on 10/11/01

% Compiler flag:
%#realonly

% Some checks
if sum(size(x1)~=size(x2)) % different sizes?
    error('NORMPDF2: The sizes of X1 and X2 must be identical!');
end;

[m n]=size(mu);

if ((m~=1)|(n~=2))&((m~=2)|(n~=1)) % wrong dimensions?
    error('NORMPDF2: MU has wrong dimensions!');
end;

if m==1 % wrong orientation?
    mu=mu'; % transpose it
end;

[m n]=size(cm);

if (m~=2)|(n~=2) % wrong dimensions?
    error('NORMPDF2: CM has wrong dimensions!');
end;

if det(cm)==0 % singular matrix?
    error('NORMPDF2: The covariance matrix must not be singular!');
end;

if (det(cm)<0)|(cm(1,1)<0)|(cm(2,2)<0)|(cm(1,2)~=cm(2,1)) % invalid matrix?
    error('NORMPDF2: Invalid covariance matrix!');
end;

% Calculation
y=[];
f=.5/pi/sqrt(det(cm));
cm_=cm^-1;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        x=[x1(i,j) x2(i,j)]';
        y(i,j)=f*exp(-.5*(x-mu)'*cm_*(x-mu));
    end;
end;

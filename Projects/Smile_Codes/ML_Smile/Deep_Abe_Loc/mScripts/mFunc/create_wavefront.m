function waveFront = create_wavefront(p, coeffs,r, theta, nFlag, conType)
% CREATWAVEFRONT is to creat wavefront based on Zernike coefficients

if(nargin<=5)
    conType = 'ANSI'; % defult as OSA convention
end

if length(p)~=length(coeffs)
    error('creatwavefront:NMlength','p and coeffs must be the same length.')
end

[n, m] = zernidx2nm(p, conType);

switch nargin
    case 4
        z = zernfun(n,m,r,theta);
    case 5
        z = zernfun(n,m,r,theta,nFlag);
    case 6
        z = zernfun(n,m,r,theta,nFlag);
    otherwise
        error('zernfun2:nargin','Incorrect number of inputs.')
end

waveFront = zeros(size(r),'single');
for i = 1:length(p)
    waveFront = waveFront + coeffs(i)*z(:,i);
end

    

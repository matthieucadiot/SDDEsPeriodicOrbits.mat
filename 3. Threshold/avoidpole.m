function check = avoidpole(p, imageboundary)
% Half-plane condition (118) of Lemma 23: check = true if the intervals imageboundary, covering
% the image h(dS_{bar rho}) of the boundary of the strip, lie in an open half-plane
% {y : Re(e^{-i phi}(y - p)) > 0}, so that h does not take the value p on S_{bar rho}.

if any(isnan(imageboundary))
    check = false;
else
    imageboundaryshifted = imageboundary - p; % p moved to 0
    disttopole = min(abs(imageboundaryshifted));
    if disttopole > 0
        % rotation putting the midpoint of the first interval on the positive real axis
        phi = intval(angle(mid(imageboundaryshifted(1))));
        imageboundaryrotated = exp(-1i*phi) * imageboundaryshifted;
        % arguments of all the points (angle is not available for intervals)
        angles = imag(log(imageboundaryrotated));
        % all the arguments in an interval of length < pi: (118) holds
        check = max(angles) - min(angles) < intval('pi');
    else
        check = false;
    end
end

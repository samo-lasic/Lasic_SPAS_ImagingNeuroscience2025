
function out = extract_tuned_LTE_3D(f, PS, trPS, ODF, g, q, dt, N0, parameters)
% Find the tuned projection on the grid defined by vectors in the ODF.

u = [ODF.x ODF.y ODF.z];

% limit frq range
ind = find(trPS<parameters.thresh);
f = f(ind);
PS = PS(ind,:,:);

% interpolate to speed up
[f, PS] = interpolate_power_spectra(f,PS, 1000);
tracePS = PS(:,1,1)+PS(:,2,2)+PS(:,3,3);

PSu = real(tensor_projection(PS,u));

% normalize to trace
trace_norm = sum(tracePS);
PSu = PSu/trace_norm; % mean(sum(PSu,2)) = 1/3
tracePS = tracePS/trace_norm;

% tracePSu = mean(PSu);
% figure(1),clf, hold on, plot(3*tracePSu,'-r'), plot(tracePS,'--b')

D0 = parameters.contour.D0;

if parameters.contour.restricted
    R = (D0^2/parameters.contour.D2Rm4)^(1/4);
    Dw = DwSpherical(2*pi*f,R,D0,0,50)'; % spherical restriction
else
    % FT of the exponential autocorrelation exp(-|t|/tau) -> 2a/(a^2+w^2), a = 1/tau
    tau = (1/parameters.contour.D2Rm4)^(1/2);%tau = 0.001;
    a = 1./tau;
    Dw = 2*a./(a^2+(2*pi*f').^2);
    Dw = D0*Dw/max(Dw);
end

% figure(1),clf, hold on, plot(tracePS/sum(tracePS)), plot(Dw/sum(Dw))

trD = sum(tracePS.*Dw)/D0;
Du = sum(3*PSu'.*repmat(Dw,1,size(PSu,1)))/D0;

out.tuning.v = abs(Du-trD)/trD;
out.tuning.m = mean(out.tuning.v);
out.tuning.std = std(out.tuning.v);
out.tuning.rstd = out.tuning.std/out.tuning.m;
%tuning = (tuning-min(tuning))/(max(tuning)-min(tuning));
%tuning = map_matrix_to_range(abs(Du-trD),0,1);
tuning = map_matrix_to_range(out.tuning.v,0,1);

%----  tuning index (just tuning, no other parameters)
[~,tuning_ind] = min(tuning);
tuning_ind = tuning_ind(1);
out.tunedLTE.tuning = out.tuning.v(tuning_ind);
out.tunedLTE.vec = u(tuning_ind,:);
out.tunedLTE.g = sqrt(3)*g*out.tunedLTE.vec';
out.tunedLTE.q = sqrt(3)*q*out.tunedLTE.vec';
out.tunedLTE.maxG = max(abs(out.tunedLTE.g));
out.tunedLTE.maxSlew = max(abs(diff(out.tunedLTE.g)))/dt;

%-------  tuning indices below the first contour line -------------

inds = find(tuning<parameters.contour.limits(1));

u = u(inds,:);
out.tuning_cont.amplitude1 = get_magnitude(g,u,dt);
out.tuning_cont.slew_rate1 = get_slew_rate(g,u,dt);

amplitude = map_matrix_to_range(out.tuning_cont.amplitude1.v,0,1);
slew_rate = map_matrix_to_range(out.tuning_cont.slew_rate1.v,0,1);


out.opt = 'amp x slew';
[~,ind] = min(amplitude.*slew_rate);
ind = ind(1);
out.optLTE.vec = u(ind,:);

PSu = real(tensor_projection(PS,out.optLTE.vec));
% normalize to trace
PSu = PSu/trace_norm; % mean(sum(PSu,2)) = 1/3
Du = sum(3*PSu'.*repmat(Dw,1,size(PSu,1)))/D0;
out.optLTE.tuning = abs(Du-trD)/trD;

out.optLTE.g = sqrt(3)*g*out.optLTE.vec';
out.optLTE.q = sqrt(3)*q*out.optLTE.vec';

out.optLTE.maxG = max(abs(out.optLTE.g));
out.optLTE.maxSlew = max(abs(diff(out.optLTE.g)))/dt;

end


function vecRotMatDemo()
numSamples = 1000000;
%% Test One:
tic;
fprintf(1,'%s',[mfilename,' :: Conducting Random Vector Rotation Test ... ']);
vec1 = randn(numSamples,3);
vec1 = bsxfun(@rdivide,vec1,sqrt(sum(vec1.^2,2)));
vec2 = randn(numSamples,3);
vec2 = bsxfun(@rdivide,vec2,sqrt(sum(vec2.^2,2)));
R = vecRotMat(vec1,vec2);
fprintf(1,'%s\n',[num2str(toc),'(sec)']);
%Test to ensure that R*vec1' = vec2'
sol = computeSol(R,vec1);
%Compute and Report Maximum Error in Rotation Matrices   
checkResults(sol,vec2);
%% Test Two:
tic;
fprintf(1,'%s',[mfilename,' :: Conducting Closely Parallel Rotation Test ... ']);
vec1 = randn(numSamples,3);
vec1 = bsxfun(@rdivide,vec1,sqrt(sum(vec1.^2,2)));
vec2 = vec1 + eps.*round(randn(numSamples,3)*10);
vec2 = bsxfun(@rdivide,vec2,sqrt(sum(vec2.^2,2)));
R = vecRotMat(vec1,vec2);
fprintf(1,'%s\n',[num2str(toc),'(sec)']);
%Test to ensure that R*vec1' = vec2'
sol = computeSol(R,vec1);
%Compute and Report Maximum Error in Rotation Matrices   
checkResults(sol,vec2);
end
function sol = computeSol(R,vec1)
sol = [squeeze(sum(bsxfun(@times,R(1,:,:),permute(vec1,[3 2 1])),2)), ...
       squeeze(sum(bsxfun(@times,R(2,:,:),permute(vec1,[3 2 1])),2)), ...
       squeeze(sum(bsxfun(@times,R(3,:,:),permute(vec1,[3 2 1])),2))];
end
function checkResults(sol,vec2)
RMS = sqrt(sum((sol-vec2).^2,2)); 
disp(['MAX RMS ERROR:' num2str(max(RMS))]);
err = abs(sol-vec2);
if max(err(:)) <= 1e-13
   disp(['TEST PASSED',' Maximum Error is ',num2str(max(err(:)))]);
else
   error(['TEST FAILED',' Maximum Error is ',num2str(max(err(:)))]);
end
end

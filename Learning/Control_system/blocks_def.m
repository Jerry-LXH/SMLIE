numg=[1]; deng=[500 0 0]; sysg=tf(numg, deng);
numh=[1 1]; denh=[1 2]; sysh=tf(numh, denh);
sys=series(sysg,sysh); 
sys=parallel(sysg,sysh);
sys=feedback(sysg,sysh,-1); % Negative feedback
sys
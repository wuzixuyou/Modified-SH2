function ML = ML_Phase(ML)

Lenses=cell(1,length(ML.coor));
for j=1:length(Lenses)
    S = ML.Pupil.*exp(-1i*SH.k.*Lenses{j}); %complex phase screen
end
ML.phase = S;
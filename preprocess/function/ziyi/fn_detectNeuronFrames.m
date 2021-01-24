function [nNeuron, nFrames, TC] = fn_detectNeuronFrames(TC)
if size(TC,1) >= size(TC,2)
    nNeuron =  size(TC,2); nFrames = size(TC,1);
    TC = TC'; 
    disp('Data Originally Arranged in Frames X Neuron')
    disp('Reshape Data into Neuron X Frames')
else
    nNeuron =  size(TC,1); nFrames = size(TC,2);
end


end
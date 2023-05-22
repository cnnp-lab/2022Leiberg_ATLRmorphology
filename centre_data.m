function [dataPre, dataPost] = centre_data(dataPre, dataPost)
% Centre data subject-wise

    dataMean = (dataPre + dataPost)/2;
    dataPre = dataPre - dataMean;
    dataPost = dataPost - dataMean;

end
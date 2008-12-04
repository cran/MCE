`calcPB` <-
function(data, index, truth) (mean(data[index]) - truth) / truth * 100


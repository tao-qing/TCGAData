path="/home/tao/taoqing/projects/sgi/TCGAArray6CNVEstimate/"
cat tumorsample_information.txt | while read line
 do
     oldname=$(echo $line | awk '{print $1"_cnvGeneMat.txt"}');
     newname=$(echo $line | awk '{print $2"_"$3"_"$4"_cnvGeneMat.txt"}'); 
    cp ${path}${oldname}  $PWD/${newname}

 done

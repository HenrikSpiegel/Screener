cat combined_bgc.fa | awk '{
        if (substr($0, 1, 1)==">") {
            st=index($0," ")
            filename=(substr($0,2,st-2) ".fasta") #takes only the first "word" after > ie ">acc.version some more description" -> "acc.version"
            }
        print $0 >> filename
        close(filename)
}'

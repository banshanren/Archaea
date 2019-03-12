library(pracma)

CalculatHurst<-function(inPath,outPath) {
    lines<-readLines(inPath)
    flag<-1
    for(line in lines){
        #convert a dna seq like '12334' which is string into a numeric: c(1 2 3 3 4)
        line<-strsplit(line,split='')                     #R auto ignore '\n', thus doesn't need strip()
        line<-as.numeric(unlist(line))

        #compute hurst exponent
        result<-hurstexp(line,display=FALSE)
        if (flag==1){
            write.table(result, file = outPath, row.names = F, quote = F, sep="\t") 
        }
        else{
            write.table(result, file = outPath, row.names = F, col.names = FALSE, quote = F, sep="\t", append = TRUE) 
        }
        flag<-flag+1
    }
}

CalculatHurst('./data/hurst/archaea_pos80_num1.txt','./data/hurst/hurst_pos_1.txt')
CalculatHurst('./data/hurst/archaea_pos80_num2.txt','./data/hurst/hurst_pos_2.txt')
CalculatHurst('./data/hurst/archaea_pos80_num3.txt','./data/hurst/hurst_pos_3.txt')
CalculatHurst('./data/hurst/archaea_neg80_num1.txt','./data/hurst/hurst_neg_1.txt')
CalculatHurst('./data/hurst/archaea_neg80_num2.txt','./data/hurst/hurst_neg_2.txt')
CalculatHurst('./data/hurst/archaea_neg80_num3.txt','./data/hurst/hurst_neg_3.txt')

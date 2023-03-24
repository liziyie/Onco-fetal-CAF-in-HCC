# Color Panel----
library(ggsci)

c102 <- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941", #1
          "#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6", #2
          "#63FFAC","#B79762","#004D43","#8FB0FF","#997D87", #3
          "#5A0007","#809693","#6A3A4C","#1B4400","#4FC601", #4
          "#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900", #5
          "#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA", #6
          "#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299", #7
          "#300018","#0AA6D8","#013349","#00846F","#372101", #8
          "#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2", #9
          "#C2FF99","#001E09","#00489C","#6F0062","#0CBD66", #10
          "#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66", #11
          "#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459", #12
          "#456648","#0086ED","#886F4C","#34362D","#B4A8BD", #13
          "#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F", #14
          "#938A81","#575329","#00FECF","#B05B6F","#8CD0FF", #15
          "#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7", #16
          "#A77500","#6367A9","#A05837","#6B002C","#772600", #17
          "#D790FF","#9B9700","#549E79","#FFF69F","#201625", #18
          "#72418F","#BC23FF","#99ADC0","#3A2465","#922329", #19
          "#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98", #20
          "#A4E804","#324E72")                               #21

c28 <- c("#023fa5","#7d87b9","#bec1d4","#d6bcc0","#bb7784", #1
         "#8e063b","#4a6fe3","#8595e1","#b5bbe3","#e6afb9", #2
         "#e07b91","#d33f6a","#11c638","#8dd593","#c6dec7", #3
         "#ead3c6","#f0b98d","#ef9708","#0fcfc0","#9cded6", #4
         "#d5eae7","#f3e1eb","#f6c4e1","#f79cd4","#7f7f7f", #5
         "#c7c7c7","#1CE6FF","#336600")                     #6 

c54 <- c("#C9BDB2","#EBDBE4","#63B472","#F6C985","#F7DDD4", #1
         "#DEEAB1","#588198","#D5E7F7","#BB4A94","#ECAFCF", #2
         "#F8BFAF","#2D563D","#CAA57D","#9E6BAB","#AFB2B7", #3
         "#6B6A6B","#66CEF6","#74517B","#F3746C","#D6D4EB", #4
         "#EACF68","#F8F4A8","#B9A96B","#EF5276","#EAA944", #5
         "#69B4CE","#7A8DBB","#A0D7C9","#C35338","#86382A", #6
         "#E3CFE6","#CBEBF6","#497C76","#7673AE","#BDB1A3", #7
         "#CFE3C6","#7AAF93","#A79388","#47885E","#409ECC", #8
         "#EAA86C","#B67464","#3955A1","#874471","#A97443", #9
         "#D69782","#D14C4C","#6A5648","#9F82A4","#97D174", #10
         "#C13B5F","#E18772","#CBB7DF","#75CAB5")           #11

c7 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

c_pat6 <- c("#e8dae3","#de5178","#e4a9d1","#9c5cad","#72477d",
            "#a9d8c7","#6ebb6e","#daf2ae","#33583a")

c_ndcirco <- c("#fb942a","#d945a0","#415f96","#c8ebe8","#f2a7b8",
               "#516dd3","#f5f179","#fbc160","#5cbac8","#97c6df",
               "#f97d3c","#8457c0","#dbc5d5","#7437c7","#e93e83")

c_odditynaive <- c("#cb8afe","#f6a9a1","#eec0da","#fd917a","#6b23e3","#ffc8c9")

Binomial_color_panel <- c("TRUE" = "#E64B35", "FALSE" = "lightgrey")

Tissue_color_panel <- c("Tumor" = "#EC7000", "Adj Normal" = "#019FAA", "Normal" = "#019FAA", "Fetal" = '#F0CE39', "Healthy" = "#9E6BAB") 

# Easy seperated color panel----
library(Polychrome)
CreateColorPanel <- function(label, rand.seed = 903417, color.seed = c("#4DCDF0", "#32CA32"), M = 1000){
  set.seed(rand.seed)
  color_panel <- createPalette(length(unique(label)), color.seed, M)
  names(color_panel) <- unique(label)
  return(color_panel)
}
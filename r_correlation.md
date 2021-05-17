# THIS IS FOR GETTING PLUS COLLELATION VALUS OF THE DATA.TXT FILE

```
library(psych)

myData <- read.table("C:/Users/Snijesh/Desktop/input.txt", header=T, sep="\t")

myData1 <- myData[-1]
rownames(myData1) <- myData[,1]
Corrt <- cor(t(myData1))

Corrt[Corrt >= 0.5]
is.na(Corrt) <- Corrt < 0.5
Corrt

#FOR GETTING THE PAIRED ITEMS
zdf <- as.data.frame(as.table(Corrt))

#FOR SAVING THE PAIRED CORRELATION VALUES
write.table(zdf, file='C:/Users/Snijesh/Desktop/output__corr_plus.txt', quote=FALSE)
```
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# THIS IS FOR GETTING MINUS COLLELATION VALUS OF OF THE DATA.TXT FILE

```
library(psych)
myData <- read.table("C:/Users/Snijesh/Desktop/input.txt", header=T, sep="\t")

myData1 <- myData[-1]
rownames(myData1) <- myData[,1]
Corrt <- cor(t(myData1))

Corrt[Corrt <= -0.5]
is.na(Corrt) <- Corrt > -0.5
Corrt

#FOR GETTING THE PAIRED ITEMS
zdf <- as.data.frame(as.table(Corrt))


#FOR SAVING THE PAIRED CORRELATION VALUES
write.table(zdf, file='C:/Users/Snijesh/Desktop/output_corr_minus.txt', quote=FALSE)
```

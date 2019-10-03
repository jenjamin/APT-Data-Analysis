# Creating mock pos file 400000 ions 
# 20 x 20 x 25 nm

zLength = 25
HalfzLength = zLength/2

x <- sample(seq(0, 20, 0.00001), 400000, replace = FALSE)
y <- sample(seq(0, 20, 0.00001), 400000, replace = FALSE)
z <- sample(seq(-HalfzLength, HalfzLength, 0.00001), 400000, replace = FALSE)

FeMatrix = 98.99
FeGB = 86
FeDiff = FeMatrix - FeGB

NiMatrix = 1
NiGB = 8
NiDiff = NiMatrix - NiGB

PMatrix = 0.01
PGB = 6
PDiff = PMatrix - PGB

SimulatedPos <- data.frame(x,y,z) %>%
  mutate(m = if_else(
    z < -1 | z > 1,
    sample(c(28,29,15.5), n(), prob = c(FeMatrix,NiMatrix,PMatrix), replace = TRUE),
    sample(c(28,29,15.5), n(), prob = c(rnorm(1,FeMatrix-(FeDiff-FeDiff*abs((z)/HalfzLength)),1),
                                      rnorm(1,NiMatrix-(NiDiff-NiDiff*abs((z)/HalfzLength)),1),
                                      rnorm(1,PMatrix-(PDiff-PDiff*abs((z)/HalfzLength)),1)),
           replace = TRUE))
  )

IonList <- lapply(SimulatedPos$m, function(MassToCharge) RangesDF2$Ion[between(MassToCharge,
                                                                               RangesDF2$Start,
                                                                               RangesDF2$End)])
IonList <- lapply(IonList, function(x) if(identical(x, character(0))) NA_character_ else x)
SimulatedPos$Ion <- unlist(IonList)
rm(IonList)

ggplot(SimulatedPos %>%
         filter(m == 15.5),
       aes(x = z)) +
  stat_bin(aes(y=cumsum(..count..)),geom="step") +
  xlim(-12.5,12.5)

a <- SimulatedPos %>% 
  group_by(Distance = cut(z, breaks= seq(-HalfzLength, HalfzLength, by = 0.1)),
           Ion) %>%
  summarise(Ioncount = n()) %>%
  ungroup() %>%
  spread(Ion, Ioncount) %>%
  mutate(Distance = as.numeric(as.character(zLength*(row_number()/n()))),
         P = replace_na(P, 0))

ggplot(a) +
  geom_point(aes(Distance, Fe))

c(FeMatrix-(FeDiff-FeDiff*abs((z)/HalfzLength)) +
    NiMatrix-(NiDiff-NiDiff*abs((z)/HalfzLength)) +
    PMatrix-(PDiff-PDiff*abs((z)/HalfzLength)))

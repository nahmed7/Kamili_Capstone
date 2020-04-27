```
library(tidyverse)
library(pwr)
library(treemapify)
library(ggformula)
library(lsr)
library(DescTools)
```
## 1. Significance: 

Many microbes evade host immunity by expressing host antigens on the surface, a phenomenon known as molecular mimicry. For example, Groub B Streptococcus (GBS), the leading cause of neonatal meningitis and sepsis, expresses a heavily sialylated outer capsule. Sialic acid is an immune-inhibitory antigen abundantly expressed on mammalian cells that serves to inhibit immune activation upon binding, thereby protecting from autoimmunity. GBS is known to exploit this mechanism to render the mammalian immune system inert. Recent studies indicate that while sialylation-mediated immune inhibition by GBS holds true for adaptive immune cells, an innate immune lectin known as Galectin-8 (Gal-8) compromises GBS viability in vitro. However, it is not known whether Gal-8 likewise plays a role physiologically in the immune response against GBS. 
**GBS is a leading cause of meningitis, pneumonia, sepsis and death in neonates. One in four women are GBS carriers, risking transmission to the baby in utero or during birth. Given in vitro findings of Gal-8 activity against GBS, we ask "does Gal-8 provide immunity in murine FRT colonized with GBS?"**

## 2. Unknown: 

It is unknown whether **Galectin-8 provides immunity in the female reproductive tract (FRT) colonized with GBS.**

## 3. Prediction:

To determine whether Gal-8 provides immunity against GBS in the murine FRT, we will measure GBS burden (colony forming units per mL of sample) in wild-type and Gal-8 knockout mice 72 hours after vaginal innoculation with GBS. Measurements taken in a pilot study were used in simulations below. 
**If Gal-8 compromises GBS viability, then Gal-8 knockout mice will have a greater GBS burden.**


## 4. Variables:

**The dependent variable** will be GBS burden, calculated as colony forming units (CFU) of GBS on GBS-selective agar plates plated with FRT sample supernatant of each mouse. The dependent variable is a continous, measured variable.
**The independent variable** is the physiological presence or absence of Gal-8. The independent variable is a discrete, sorted variable of two levels, wild-type and Gal-8 knockout.
The experimental mice will be commercially sourced and age matched. They will not be from the same litter and will be treated as independent from one another.

## 5. Statistical Hypothesis:

The null hypothesis is that GBS burden in Gal-8 knock out FRT is less than or equal to the GBS burden in wild-type FRT. The alternate hypothesis is that GBS burden in Gal-8 knock out FRT is greater than the GBS burden in wild-type FRT. 

**The null hypothesis is that Gal-8 is not involved in immunity against GBS in the murine FRT.** 
**Null Hypothesis:** µ(GBS_WT) greater than or equal to µ(GBS_Gal-8KO)

**The alternate hypothesis is that Gal-8 is involved in immunity against GBS in murine FRT.**
**Alternate Hypothesis:** µ(GBS_WT) less than µ(GBS_Gal-8KO)


## 6. Statistical Test:

We will employ an **unpaired, one-sided t-test** for comparing group means of GBS burden in wild-type and knockout mice. The unpaired, one-sided t-test is appropriate in this case because the mice are treated as independent of each other and the response variable is a continuous, measured variable. Furthermore, the null hypothesis is one-sided: GBS burden in the wild-type group is greater than or equal to the Gal-8 knockout group. The sample size will be based upon a **power of 80%**, tolerace for **type2 error will be 20%**. The decision threshold for **type1 error will be 5%**. The null hypothesis will be rejected at a p-value of less than 0.05.

## 7. Procedures:

8-week old female wild-type and Gal-8 KO mice will be commercially sourced  from Jackson Labs. Each mouse will recieve beta-estradiol to synchronize their estrous cycle, as previous studies indicate that varying stages of estrous contributes to variations in infection levels. Mice will then be innoculated with 10^7 GBS cells via vaginal innoculation 24 hours after beta-estradiol administration. Mice will be sacrificed and FRT will be harvested, washed and vortexed in PBS. The supernatant will be plated onto GBS selective agar plates 72 hours later. GBS burden will be measured by quantification of GBS CFU/mL. 
Data will be analyzed by an unpaired, one-sided student’s t-test to determine statistical differences between the two groups. We will use a type 1 alpha error threshold of 0.05 and aim for a statistical power of 0.80 with a resulting type 2 beta of 0.20. The number of mice per group to achieve the desired alpha and beta threshold will be calculated using a power analysis in R, as well as a Monte Carlo simulation.

## 8. Graph:

```
n <- 5
mean_WT <- 50000
sd_WT <- 20000
mean_KO <- 500000
sd_KO <- 200000

data <- tibble(WT=abs(rnorm(n, mean_WT, sd_WT)),
         Gal8_KO=abs(rnorm(n, mean_KO, sd_KO))
         
         ) %>% 
    
    pivot_longer(cols=WT:Gal8_KO, 
                 names_to="Group", 
                 values_to="GBS_CFU")

data
```
<img width="306" alt="Screen Shot 2020-04-27 at 3 22 43 PM" src="https://user-images.githubusercontent.com/64428885/80412120-55433e00-889b-11ea-9634-a08b4aca9715.png">

```{r}
ggplot(data, aes(x=Group, y=GBS_CFU))+
    stat_summary(
    fun.data=mean_sdl,
    fun.args = list(mult=1),
    geom="crossbar",
    width=0.4,
    color="red"
  )+
  geom_jitter(aes(Group, GBS_CFU))+
  theme_classic()+
  labs(y="GBS Burden (CFU/mL)", title="GBS Burden in FRT of Wild-type and Galectin-8 KO Mice")

```
![Task8plot.pdf](https://github.com/nahmed7/Kamili_Capstone/files/4541631/Task8plot.pdf)

## 9. Sample Size:

**Based on an alpha of 0.05 and beta of 0.20, we will need a minimum of 4 mice per group for a sufficiently powered experiment per the power analysis and the monte carlo simulation below:**

```
t.pwr <- function(n){
  # Intitializers
  mWT=50000; sdWT=20000; mG8KO= 500000; sdG8KO=200000
  
  ssims=1000
  p.values <- c()
  i <- 1
  
  #THE monte carlo
  repeat{
    x=rnorm(n, mWT, sdWT); 
    y=rnorm(n, mG8KO, sdG8KO);
    p <- t.test(x, y, 
                paired=F, 
                alternative="less", 
                var.equal=F,
                conf.level=0.95)$p.value
    p.values[i] <- p
    if (i==ssims) break
    i = i+1
    pwr <- length(which(p.values<0.05))/ssims
  }
  return(pwr)
}


frame <- data.frame(n=2:50)
data <- bind_cols(frame, 
                  power=apply(frame, 1, t.pwr))
#plot
ggplot(data, aes(n, power))+
  geom_point() +
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  scale_x_continuous(breaks=c(seq(0,50,2)))+
  labs(x="n per group")

```

```

alpha <- 0.05
nSims <- 10000
p <-numeric(nSims) 

# the monte carlo function

for(i in 1:nSims){ 
  WT <-rnorm(n=4, mean=500000, sd=200000)
  Gal8KO <-rnorm(n=4, mean=50000, sd=20000)
  z<-t.test(WT,Gal8KO, paired=FALSE, alternative="greater") 
  p[i]<-z$p.value 
}

hits <- length(which(p < alpha));hits
power <- hits/nSims;power

#main="Histogram of p-values under the null",
hist(p, col = "blue", ylim = c(0, 1000), xlim = c(0.0, 1.0), main ="Histogram of simulated p-values", xlab=("Observed p-value"))

```

#Anogenital warts (AW) incidence per age and gender

#From RIVM raport sexually transmitted infections 2021 Table 7.1
#We set the averages below and above 25 to these numbers:
#women <25: 2.7, women >25: 2.1, men <25: 2.1, men >25: 3.5

#From Woestenberg et al. CID 2020 
#Incidence for women age 12-22.
GW_inc1222_w <- c(0.14,0.09,0.38,0.57,1.57,3.19,4.24,5.57,5.19,8.19,6.80)
#We assume a linear trend from age 22 to 25

#From Seksuele gezondheid in Nederland 2011 Table 6: Trend sexual partners
#women
#25-39 (5.9): x
#40-54 (4.4): 0.746x
#55-70 (2.0): 0.339x
#men
#25-39 (13.4): x
#40-54 (11.7): 0.873x
#55-70 (6.4): 0.478x

#From leefstijlmonitor seksuele gezondheid 2019 Table 2: trend 70-85 
#55-64: 3.4
#75-85: 0.4
#women 70-85: 0.339*0.118 = 0.04 
#men 70-85: 0.478*0.118 = 0.056

trend_vec_w <- c(1,0.746,0.339,0.04)
trend_vec_m <- c(1,0.873,0.478,0.056)

#We assume 0 incidence before age 12 and after age 84,
#and take this into account when computing the averages

#Let x be the incidence at 25
AverageAbove25 <- function(x,trend_vec){
  x*sum(trend_vec*c(15,15,15,15))/76
}

f <- function(x,trend_vec,av.val){
  AverageAbove25(x,trend_vec) - av.val
}

#We want the average for women above 25 equal to 2.1
sol.x_w <- uniroot(f,interval=c(0,20),trend_vec_w,2.1)
Inc25_w <- sol.x_w$root
Inc_w <- trend_vec_w*Inc25_w

#We want the average for men above 25 equal to 3.5
sol.x_m <- uniroot(f,interval=c(0,20),trend_vec_m,3.5)
Inc25_m <- sol.x_m$root
Inc_m <- trend_vec_m*Inc25_m

#Let y be the incidence at 12
TrendUnder25 <- GW_inc1222_w/GW_inc1222_w[1]

AverageUnder25 <- function(y,Inc25){
  Inc22 <- TrendUnder25[11]*y
  (sum(TrendUnder25*y)+ (Inc22-(Inc22-Inc25)/3)+(Inc22-2*(Inc22-Inc25)/3) )/24
}

g <- function(y,inc25,av.val){
  AverageUnder25(y,inc25) - av.val
}

#We want the average for women below 25 equal to 2.7
sol.y <- uniroot(g,interval=c(0,20),Inc25_w,2.7)
Inc12_w <- sol.y$root
Inc1222_w <- TrendUnder25*Inc12_w
Inc23_w <- Inc1222_w[11]-(Inc1222_w[11]-Inc25_w)/3
Inc24_w <- Inc1222_w[11]-2*(Inc1222_w[11]-Inc25_w)/3

#We want the average for men below 25 equal to 2.1
sol.y <- uniroot(g,interval=c(0,20),Inc25_m,2.1)
Inc12_m <- sol.y$root
Inc1222_m <- TrendUnder25*Inc12_m
Inc23_m <- Inc1222_m[11]-(Inc1222_m[11]-Inc25_m)/3
Inc24_m <- Inc1222_m[11]-2*(Inc1222_m[11]-Inc25_m)/3

#Incidence per gender for ages 1:100
AW.inc_w <- c(numeric(11),Inc1222_w,Inc23_w,Inc24_w,
              rep(Inc_w,c(15,15,15,15)),
              rep(0,16))

AW.inc_m <- c(numeric(11),Inc1222_m,Inc23_m,Inc24_m,
              rep(Inc_m,c(15,15,15,15)),
              rep(0,16))

#plots
plot(1:100,AW.inc_m,type="l",col="#00B6EB",main="AW episodes per 1000 p-yrs",ylim=c(0,12),
     ylab="",xlab="age")
lines(AW.inc_w,type="l",col="#F8766D")
legend('topright', legend=c("Women", "Men"),
       col=c("#F8766D", "#00B6EB"), lty=1, cex=.9,
       text.font=4)

#Check:
mean(AW.inc_w[1:24]) #2.7
mean(AW.inc_w[25:100]) #2.1

mean(AW.inc_m[1:24]) #2.1
mean(AW.inc_m[25:100]) #3.5  

library(readr)
# test case 1 marstestdata.rda
load("~/MARS/MARS/data/marstestdata.rda")
y <- marstestdata$y
fit1 <- lm(y~.-1, marstestdata)
out1 <- mars(fit1,marstestdata,control=mars.control())
print.mars(out1)
summary.mars(out1)
predict.mars(out1)
anova.mars(out1)
plot.mars(out1)

# test case 2 laptop prices
laptop_price <- read_csv("data/laptop_price.csv")
y <- laptop_price$Price_euros
fit2 <- lm(y~.-1, laptop_price)
out2 <- mars(fit2,laptop_price,control=mars.control())
print.mars(out2)
summary.mars(out2)
predict.mars(out2)
anova.mars(out2)
plot.mars(out2)

# test case 3 real estate data
Real_estate <- read_csv("data/Real estate.csv")
y <- Real_estate$`Y house price of unit area`
fit3 <- lm(y~.-1, Real_estate)
out3 <- mars(fit3,Real_estate,control=mars.control())
print.mars(out3)
summary.mars(out3)
predict.mars(out3)
anova.mars(out3)
plot.mars(out3)

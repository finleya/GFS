h=trees$Height
d=trees$Girth
v=trees$Volume
d2h = d*d*h
par(mai=c(1,1,0.2,1))
par(mfrow=c(1,1))
plot(d2h,v,ylab = expression(paste("V"," (","ft"^"3",")")),
     xlab= expression(paste( "D"^"2","H"," (in"^"2","ft)")))


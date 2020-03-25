total.max.add = 350

png("40-30_sample_size_estimation.png")
#15% versus 5%
print("15% versus 5%")
large.start = 40
small.start = 30

total.size = c()
pvalues = c()

for (i in 0:total.max.add){
	test.large = large.start + i
	test.small = small.start + i
	
	integer.pos.large = as.integer(0.15 * test.large)
	integer.pos.small = as.integer(0.05 * test.small)
	
	#2-sided, since you don't know which is larger ahead of time
	fisher.mat = matrix(c(integer.pos.large, test.large, integer.pos.small, test.small), ncol=2)
	rownames(fisher.mat) = c("pos","neg")
	colnames(fisher.mat) = c("higher","lower")
	
	result = fisher.test(fisher.mat)
	
	total.size = c(total.size, test.large + test.small)
	pvalues = c(pvalues, result$p.value) 
}#end for (i in 0:total.max.add)

plot(total.size, pvalues,
	ylim=c(0,1), xlim=c(0,max(total.size) + 5),
	col=rgb(red=0, green=0, blue=0, alpha=0.5), pch=16)
max.not.detected = max(total.size[pvalues > 0.05])
print(max.not.detected)
abline(v=max.not.detected, col = rgb(red=0, green=0, blue=0, alpha=0.8), lty=4)
abline(h=0.05)
loess.model = loess(pvalues ~ total.size)
lines(total.size, predict(loess.model), col = "black", lwd=4)
rect(-10, -0.1, max(total.size) + 100, 0.05, col=rgb(red=0, green=1, blue=0, alpha=0.2), border=NA)
legend("topright", c("15%-vs-5%","10%-vs-5%"),
		col=c("black","blue"), lwd=4)
		
#10% versus 5%
print("10% versus 5%")
large.start = 40
small.start = 30

total.size = c()
pvalues = c()

for (i in 0:total.max.add){
	test.large = large.start + i
	test.small = small.start + i
	
	integer.pos.large = as.integer(0.10 * test.large)
	integer.pos.small = as.integer(0.05 * test.small)
	
	#2-sided, since you don't know which is larger ahead of time
	fisher.mat = matrix(c(integer.pos.large, test.large, integer.pos.small, test.small), ncol=2)
	rownames(fisher.mat) = c("pos","neg")
	colnames(fisher.mat) = c("higher","lower")
	
	result = fisher.test(fisher.mat)
	
	total.size = c(total.size, test.large + test.small)
	pvalues = c(pvalues, result$p.value) 
}#end for (i in 0:total.max.add)

points(total.size, pvalues,
	col=rgb(red=0, green=0, blue=1, alpha=0.5), pch=16)
max.not.detected = max(total.size[pvalues > 0.05])
print(max.not.detected)
abline(v=max.not.detected, col = rgb(red=0, green=0, blue=1, alpha=0.8), lty=4)
loess.model = loess(pvalues ~ total.size)
lines(total.size, predict(loess.model), col = "blue", lwd=4)
rect(-10, -0.1, max(total.size) + 100, 0.05, col=rgb(red=0, green=1, blue=0, alpha=0.2), border=NA)
dev.off()
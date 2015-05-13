from subprocess import call 
import time 
import matplotlib.pyplot as plt

# currently changing NM parameter 

overlap = []
entropy = []
runningtime = []

for x in range(0, 10): 


	#default settings
	call(["./motif", "8", "1", "500", "10", str(x)])	

	sites = open("data_set" + str(x) + "/sites.txt", "r")
	predicted_sites = open("data_set" + str(x) + "/predictedsites.txt", "r")

	count = 0
	array_actual = []
	array_predicted = []
	for line in sites: 
		array_actual.append(line)

	for otherline in predicted_sites: 
		array_predicted.append(otherline)	

	for x in range(0, 10):
		if array_actual[x] == array_predicted[x]:
			count += 1	

	overlap.append(count)

	print "number of overlapping sites: " + str(count)
	print "\n"

	time.sleep(2)


f = plt.plot(overlap)
plt.ylabel('number of overlapping sites')
f.show()

e = open("entropy.txt")
r = open("runningtime.txt")

for line in e: 
	entropy.append(line)

for line in r: 
	runningtime.append(line); 

g = plt.plot(entropy)
plt.ylabel('relative entropy')
g.show()

h = plt.plot(runningtime)
plt.ylabel('running times')
h.show()

raw_input()



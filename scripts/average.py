def average(mylist):
	n = len(mylist)
	sum = 0
	for i in mylist:
		sum += i
	return sum/n

mylist = [1,2,3,4,5,6,7,8,9,10]
print(average(mylist))
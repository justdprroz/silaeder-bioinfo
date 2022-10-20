import csv
import sys
import matplotlib.pyplot as plt

l = []

reader = csv.reader(sys.stdin, delimiter=',')
for row in reader:
    l.append(float(row[3]))

plt.hist(l, bins=100)
plt.show()
import re
import sys
import csv


def verify_timings(output_file):
	open_file = open(output_file, 'r')
	try:
		lines = open_file.readlines()
	except UnicodeDecodeError:
		file.close()
		open_file = open(folder + file, encoding = "ISO-8859-1")
		lines = open_file.readlines()
	counter = 1
	for line in reversed(lines):
		if counter > 30:
			return -1
		else:
			if "elapsed time" in line: # this is seconds
				hours = re.findall("\d+\.\d+", line)
				return hours[0]
			counter += 1

celia = ("celia", "./test_cases/output/celia1990/celia.txt")
colbeck_1 = ("colbeck_1", "./test_cases/output/colbeck1976/colbeck_1.txt")
colbeck_2 = ("colbeck_2", "./test_cases/output/colbeck1976/colbeck_2.txt")
colbeck_3 = ("colbeck_3", "./test_cases/output/colbeck1976/colbeck_3.txt")
miller_clay = ("miller_clay", "./test_cases/output/miller1998/miller_clay.txt")
miller_loam = ("miller_loam", "./test_cases/output/miller1998/miller_loam.txt")
miller_sand = ("miller_sand", "./test_cases/output/miller1998/miller_sand.txt")
mizoguchi = ("mizoguchi", "./test_cases/output/mizoguchi1990/mizoguchi.txt")
vanderborght_1 = ("vanderborght_1", "./test_cases/output/vanderborght2005/vanderborght_1.txt")
vanderborght_2 = ("vanderborght_2", "./test_cases/output/vanderborght2005/vanderborght_2.txt")
vanderborght_3 = ("vanderborght_3", "./test_cases/output/vanderborght2005/vanderborght_3.txt")
wigmosta_1 = ("wigmosta_1", "./test_cases/output/wigmosta1999/wigmosta_1.txt")
wigmosta_2 = ("wigmosta_2", "./test_cases/output/wigmosta1999/wigmosta_2.txt")


all_test_cases = [celia, colbeck_1, colbeck_2, colbeck_3, miller_clay, miller_loam, miller_sand, \
	mizoguchi, vanderborght_1, vanderborght_2, vanderborght_3, wigmosta_1, wigmosta_2]
# miller_clay and miller_sand take the most amount of time
fast_test_cases = [celia, colbeck_1, colbeck_2, colbeck_3, miller_loam, mizoguchi, \
	vanderborght_1, vanderborght_2, vanderborght_3]
	
# get all the verification data
if len(sys.argv) == 1:
	print("Need to which tests to run fast or all")
	print("Specify all tests with: all")
	print("Specify fast tests with: fast")
else:
	if sys.argv[1] == "all":
		test_cases = all_test_cases
	elif sys.argv[1] == "fast":
		test_cases = fast_test_cases
	else:
		print("Please use command line arguments fast or all")


with open("timing_info.csv", 'w') as timing_info:

	writer = csv.writer(timing_info)
	header = ["test", "time_in_seconds"]

	writer.writerow(header)

	for test in test_cases:
		seconds = verify_timings(test[1])
		print(seconds)
		writer.writerow([test[0], seconds])





import os
import time
import shutil
positive_path = '..\\audiofile\\positive'
positive_files = os.listdir(positive_path)
num_p = len(positive_files)
negetive_path = '..\\audiofile\\negetive'
negetive_files = os.listdir(negetive_path)
num_n = len(negetive_files)
loopnum=num_p+num_n
counter_files=os.listdir('.\\counter')
counter=len(counter_files)
scale = 50
print("Feature Extracting...".center(scale // 2,"-"))
start = time.perf_counter()
while counter<loopnum :
	a = "*" * round(counter / loopnum * scale)
	b = "." * round(scale - counter / loopnum *scale)
	c = counter / loopnum * 100
	dur = time.perf_counter() - start
	print("\r{:^3.0f}%[{}->{}]{:.2f}s".format(c,a,b,dur),end = "")
	counter_files=os.listdir('.\\counter')
	counter=len(counter_files)
	time.sleep(0.1)
a = "*" * scale
b = "." * round(scale - scale)
c = 1 * 100
dur = time.perf_counter() - start
print("\r{:^3.0f}%[{}->{}]{:.2f}s".format(c,a,b,dur),end = "")
print("\n"+"complete!".center(scale // 2,"-"))
shutil.rmtree('.\\counter') 
time.sleep(3)
quit()
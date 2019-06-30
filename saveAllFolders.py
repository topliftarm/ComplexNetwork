import os

path = '.'

#folders = []

# r=root, d=directories, f = files
'''
for r, d, f in os.walk(path):
    for folder in d:
		print('Processing data in'+folder)
		os.system("python plotData.py "+ folder)
'''

cwd = os.getcwd()
print(cwd)
for r, d, f in os.walk(path):
	for folder in d:
		#print("python plotData.py "+ str(folder) + " "+ str(cwd)) 
		print("python plotData.py "+ '"'+ folder+ '"'+ " "+ '"'+ cwd+'"') 
		os.system("python plotData.py "+ '"'+ folder+ '"'+ " "+ '"'+ cwd+'"') 
		
        
		
		


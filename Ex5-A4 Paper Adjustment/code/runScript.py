import os
import subprocess


def main():
	sourceImagePath = "../TestData2/"
	outputImagePath = "../transformedImg/"
	for sourceImageName in os.listdir(sourceImagePath):
		args = ['./Ex5']
		args.append(sourceImagePath + sourceImageName)
		#args.append(outputImagePath + sourceImageName.strip('jpg') + 'bmp')
		subprocess.call(args)


if __name__ == "__main__":
	main()
